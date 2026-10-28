//
// Created by Gregory Johnson on 7/8/26.
//

#include "DataStructures/OptiRefMatrix.h"
#include <Rcpp.h>

#include "Adapters/CountTableAdapter.h"
#include "DataStructures/OptiMatrix.h"
#include "MothurDependencies/OneGapPairwiseDistance.h"
#include "MothurDependencies/PairwiseDistanceCalculator.h"
#include "RNG/SubSample.h"
//
//  optirefmatrix.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/3/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//


/***********************************************************************/

OptiRefMatrix::OptiRefMatrix(const std::vector<std::vector<long long>>& close,  const std::vector<std::string>& name,
const std::vector<std::string>& singleton, const double c, const int seed) : rng(seed), seed(seed), square(false),
                                                                             numFitDists(0),
                                                                             numRefDists(0),
                                                                             numRefSingletons(0),
                                                                             numFitSingletons(0),
                                                                             numBetweenDists(0),
                                                                             numSingletons(0),
                                                                             numFitSeqs(0),
                                                                             fitPercent(0) {
    closeness = close;
    nameMap = name;
    singletons = singleton;
    cutoff = c;
}

OptiRefMatrix::OptiRefMatrix(const OptiData* matrix, const CountTableAdapter& adapter, const double fP,
                             const std::string& refWeight, const int seed):rng(seed), seed(seed) {

    numFitSingletons = 0;
    numRefSingletons = 0;
    numSingletons = 0;
    numBetweenDists = 0;
    numFitDists = 0;
    numRefDists = 0;
    numFitSeqs = 0;
    refWeightMethod = refWeight;

    fitPercent = fP / 100.0;
    if (fitPercent < 0.001) { fitPercent = 0.10; Rcpp::warning("[WARNING]: fit percentage must be between 0.001 (0.1%) and 1.0 (100%). Setting to 0.10 or 10%. \n"); } //minumum of 0.1%
    else if (fitPercent > 100.0) {  Rcpp::stop("[ERROR]: fit percentage must be between 0.0001 and 100.0\n"); }

    square = false;
    std::unordered_set<std::string> noRefNamesSet;
    ReadFiles(matrix, adapter,noRefNamesSet);
}

OptiRefMatrix::OptiRefMatrix(const OptiData* matrix, const CountTableAdapter& adapter,
    std::unordered_set<std::string> accnosRefFileNames, const int seed):rng(seed), seed(seed) {

    numFitSingletons = 0;
    numRefSingletons = 0;
    numSingletons = 0;
    numBetweenDists = 0;
    numFitDists = 0;
    numRefDists = 0;
    numFitSeqs = 0;
    refWeightMethod = "accnos";

    square = false;
    ReadFiles(matrix, adapter,accnosRefFileNames);
}

OptiRefMatrix::OptiRefMatrix(const OptiData *refMatrix, const CountTableAdapter &refAdapter,
              const float fitPercent, const double cutoff, const int seed):rng(seed), seed(seed), cutoff(cutoff),
fitPercent(fitPercent){

    numFitSingletons = 0;
    numRefSingletons = 0;
    numSingletons = 0;
    numBetweenDists = 0;
    numFitDists = 0;
    numRefDists = 0;
    numFitSeqs = 0;
    refWeightMethod = "none";

    square = false;
    std::unordered_set<std::string> noRefNamesSet;
    ReadFiles(refMatrix, refAdapter, noRefNamesSet);
}

/***********************************************************************/
//Since we are extracting a subset of the seqs some reads that may not have been singletons
OptiData* OptiRefMatrix::extractRefMatrix() {
    std::set<long long> seqs;
    for (long long i = 0; i < static_cast<long long>(isRef.size()); i++) {
        if (isRef[i]) { seqs.insert(i); }
    }

    std::vector<std::string> subsetNameMap;
    std::vector<std::string> subsetSingletons;
    std::vector< std::vector<long long> > subsetCloseness;
    std::map<long long, long long> thisNameMap;
    std::map<long long, long long> nonSingletonNameMap;
    std::vector<bool> singleton; singleton.resize(seqs.size(), true);
    int count = 0;
    for (const auto& seq : seqs)  {
        long long seqNum = seq;
        thisNameMap[seqNum] = count;
        nonSingletonNameMap[count] = seqNum;

        std::vector<long long> thisSeqsCloseSeqs = getCloseSeqs(seqNum);
        for (const long long thisSeq : thisSeqsCloseSeqs) {
            //is this seq in the std::set of unfitted?
            if (seqs.count(thisSeq) != 0) { singleton[thisNameMap[seqNum]] = false; }
        }
        count++;
    }

    int nonSingletonCount = 0;
    for (long long i = 0; i < static_cast<long long>(singleton.size()); i++) {
        if (!singleton[i]) { //if you are not a singleton
            nonSingletonNameMap[i] = nonSingletonCount;
            nonSingletonCount++;
        }else { seqs.erase(nonSingletonNameMap[i]);  subsetSingletons.push_back(getName(nonSingletonNameMap[i])); } //remove from unfitted
    }
    singleton.clear();

    subsetCloseness.resize(nonSingletonCount);
    for (const auto& seqNum : seqs) {
        std::vector<long long> thisSeqsCloseSeqs = getCloseSeqs(seqNum);
        std::vector<long long> thisSeqsCloseUnFittedSeqs;
        for (long long thisSeq : thisSeqsCloseSeqs) {
            //is this seq in the std::set of unfitted?
            if (seqs.count(thisSeq) != 0) { thisSeqsCloseUnFittedSeqs.emplace_back(nonSingletonNameMap[thisNameMap[thisSeq]]); }
        }
        if (!thisSeqsCloseUnFittedSeqs.empty()) {
            subsetCloseness[nonSingletonNameMap[thisNameMap[seqNum]]] = thisSeqsCloseUnFittedSeqs;
            subsetNameMap.push_back(getName(seqNum));
        }

    }

    for (size_t i = 0; i < isSingleRef.size(); i++) {
        if (isSingleRef[i]) {
            subsetSingletons.push_back(singletons[i]);
        }
    }

    auto* unfittedMatrix = new OptiRefMatrix(subsetCloseness, subsetNameMap, subsetSingletons, cutoff, seed);

    return unfittedMatrix;
}
/***********************************************************************/
//given std::set of names, pull out their dists and create optimatrix
OptiData* OptiRefMatrix::extractMatrixSubset(const std::unordered_set<std::string>& seqs) {
        std::unordered_set<long long> seqIndexes = getIndexes(seqs);
        return extractMatrixSubset(seqIndexes);
}
/***********************************************************************/
//given matrix indexes of seqs, pull out their dists and create optimatrix
OptiData* OptiRefMatrix::extractMatrixSubset(std::unordered_set<long long> & seqs) {
    std::vector<std::string> subsetNameMap;
    std::vector<std::string> subsetSingletons;
    std::vector< std::vector<long long> > subsetCloseness;
    std::map<long long, long long> thisNameMap;
    std::map<long long, long long> nonSingletonNameMap;
    std::vector<bool> singleton; singleton.resize(seqs.size(), true);
    int count = 0;

    for (auto it = seqs.begin(); it != seqs.end(); it++) {

        long long seqNum = *it;
        thisNameMap[seqNum] = count;
        nonSingletonNameMap[count] = seqNum;

        std::vector<long long> thisSeqsCloseSeqs = getCloseSeqs(seqNum);
        for (long long thisSeq : thisSeqsCloseSeqs) {
            //is this seq in the std::set of unfitted?
            if (seqs.count(thisSeq) != 0) { singleton[thisNameMap[seqNum]] = false; }
        }
        count++;
    }

    int nonSingletonCount = 0;
    for (long long i = 0; i < static_cast<long long>(singleton.size()); i++) {
        if (!singleton[i]) { //if you are a singleton
            nonSingletonNameMap[i] = nonSingletonCount;
            nonSingletonCount++;
        }else { seqs.erase(nonSingletonNameMap[i]);  subsetSingletons.push_back(getName(nonSingletonNameMap[i])); } //remove from unfitted
    }
    singleton.clear();

    subsetCloseness.resize(nonSingletonCount);
    for (const auto seq : seqs) {

        long long seqNum = seq;

        std::vector<long long> thisSeqsCloseSeqs = getCloseSeqs(seqNum);
        std::vector<long long> thisSeqsCloseUnFittedSeqs;
        for (long long thisSeq : thisSeqsCloseSeqs) {
            //is this seq in the std::set of unfitted?
            if (seqs.count(thisSeq) != 0) { thisSeqsCloseUnFittedSeqs.emplace_back(nonSingletonNameMap[thisNameMap[thisSeq]]); }
        }

        if (!thisSeqsCloseUnFittedSeqs.empty()) {
            subsetCloseness[nonSingletonNameMap[thisNameMap[seqNum]]] = thisSeqsCloseUnFittedSeqs;
            subsetNameMap.push_back(getName(seqNum));
        }

    }

    return new OptiRefMatrix(subsetCloseness, subsetNameMap, subsetSingletons, cutoff, seed);
}
/***********************************************************************/
std::vector<long long> OptiRefMatrix::getTranslatedBins(std::vector<std::vector<std::string> > & binNames, std::vector<std::vector<long long> > & fixedBins) {

    fixedBins.clear();

    std::map<std::string, long long> nameIndexes;
    std::set<std::string> unique;
    const auto nameMapSize = static_cast<long long>(nameMap.size());
    const auto singletonSize = static_cast<long long>(singletons.size());
    const auto closenessSize = static_cast<long long>(closeness.size());

    for (long long i = 0; i < nameMapSize; i++) { //std::vector of std::string representing the sequences in the matrix from the name file.
        std::vector<std::string> thisSeqsReps; Utils::splitAtComma(nameMap[i], thisSeqsReps); //split redundant names
        if (i < closenessSize) {  nameIndexes[thisSeqsReps[0]] = i;  } //this is a sequence with distances in the matrix
        if (thisSeqsReps.size() == 1) { //you are unique
            unique.insert(thisSeqsReps[0]);
        }
    }

    for (long long i = 0; i < singletonSize; i++) {
        if (isSingleRef[i]) {
            std::vector<std::string> thisSeqsReps; Utils::splitAtComma(singletons[i], thisSeqsReps); //split redundant names
            nameIndexes[thisSeqsReps[0]] = -1;
            if (thisSeqsReps.size() == 1) { unique.insert(thisSeqsReps[0]); }
        }
    }

    for (auto & binName : binNames) { //for each OTU
        std::vector<long long> thisBinsSeqs;
        for (const auto & name : binName) { //for each sequence
            // const auto it = nameIndexes.find(name);
            if (nameIndexes.find(name) != nameIndexes.end()) {
                thisBinsSeqs.push_back(nameIndexes[name]);
            }
            // if (it == nameIndexes.end()) { }//not in distance matrix, but needs a value in fixedBins. 2 reasons for making it here: you are a redundant name in the listfile, you do not have any distances
            // else { thisBinsSeqs.push_back(it->second);  } //"name" of sequence in matrix
        }
        fixedBins.push_back(thisBinsSeqs);
    }

    return (getFitSeqs());

}
/***********************************************************************/
//assumes that i is a fitSeq
bool OptiRefMatrix::isCloseFit(const long long i, const long long toFind, bool& isFit) const {
    if (i < 0) return false;
    if (i > static_cast<long long>(closeness.size())) { Rcpp::stop("[ERROR]: index is not valid.\n");}

    bool found = false;
    if (!isRef[toFind]) { //are you a fit seq
        found = std::binary_search(closeness[i].begin(), closeness[i].end(), toFind);
        // if (closeness[i].find(toFind) != closeness[i].end()) {  //are you close
        //     found = true;
        // }
        isFit = true;
    } else { isFit = false;  }
    return found;
}
/***********************************************************************/
//does not include singletons, only reads in closeness
std::vector<long long> OptiRefMatrix::getRefSeqs() {
    const auto isRefSize = static_cast<long long>(isRef.size());

    std::vector<long long> refSeqsIndexes;
    refSeqsIndexes.reserve(isRefSize);

    for (long long i = 0; i < isRefSize; i++) {
        if (isRef[i]) { refSeqsIndexes.emplace_back(i); }
    }
    return refSeqsIndexes;
}
/***********************************************************************/
std::vector<std::string> OptiRefMatrix::getRefSingletonNames() {
    const auto isSingleRefSize = static_cast<long long>(isSingleRef.size());
    std::vector<std::string> refSeqsNames;
    refSeqsNames.reserve(isSingleRefSize);

    for (long long i = 0; i < isSingleRefSize; i++) {
        if (isSingleRef[i]) { refSeqsNames.emplace_back(singletons[i]); }
    }

    return refSeqsNames;
}
/***********************************************************************/
std::vector<long long> OptiRefMatrix::getFitSeqs() {
    const auto isRefSize = static_cast<long long>(isRef.size());

    std::vector<long long> fitSeqsIndexes;
    fitSeqsIndexes.reserve(isRefSize);

    for (long long i = 0; i < isRefSize; i++) {
        if (!isRef[i]) { fitSeqsIndexes.emplace_back(i);  }
    }
    return fitSeqsIndexes;
}

void OptiRefMatrix::UpdateSingletons(const FastaDatabase &refFastaSequences, const FastaDatabase &fitFastaSequences,
    std::vector<bool> &singletons, std::unordered_map<std::string, long long> &nameMap) {
    PairwiseDistanceCalculator* calculator = new OneGapPairwiseDistance();
    const std::vector<FastaData>& refData = refFastaSequences.GetFastaDataBase();
    const std::vector<FastaData>& fitData = fitFastaSequences.GetFastaDataBase();

    for (const auto &[refName, refSequence] : refData) {
        const long long refIndex = nameMap[refName];

        for (const auto &[fitName, fitSequence] : fitData) {
            const long long fitIndex = nameMap[fitName];
            // if (!singletons[fitIndex]) continue; // We only want to check if singletons are no considering non-singletons
            if (const double result = calculator->Execute(refSequence, fitSequence);
                result > cutoff) continue;
            closeness[refIndex].emplace_back(fitIndex);
            closeness[fitIndex].emplace_back(refIndex);
            singletons[refIndex] = false;
            singletons[fitIndex] = false;
        }
    }
    for (auto& close : closeness) {
        if (close.empty()) continue;
        std::sort(close.begin(), close.end());
    }
}

/***********************************************************************/
long long OptiRefMatrix::getNumFitTrueSingletons() {
    return std::count(isSingleRef.begin(), isSingleRef.end(), false);
}
/***********************************************************************/
long long OptiRefMatrix::getNumFitClose(const long long index) {
        long long numClose = 0;
        if (index < 0) return numClose;
        if (index > static_cast<long long>(closeness.size())) Rcpp::stop("[ERROR]: index is not valid.\n");
            //reference seqs all have indexes less than refEnd
        for (const long long it : closeness[index]) {
            if (!isRef[it]) {  numClose++; } //you are a fit seq
        }


        return numClose;
}
/***********************************************************************/
long long OptiRefMatrix::getNumRefClose(const long long index) {
        long long numClose = 0;

        if (index < 0) return numClose;
        if (index > static_cast<long long>(closeness.size())) Rcpp::stop("[ERROR]: index is not valid.\n");
            //reference seqs all have indexes less than refEnd
        for (const long long it : closeness[index]) {
            if (isRef[it]) {  numClose++; } //you are a ref seq
        }

        return numClose;
}
/***********************************************************************/
std::set<long long> OptiRefMatrix::getCloseFitSeqs(const long long index){
    std::set<long long> closeSeqs;

    if (index < 0) return closeSeqs;
    if (index > static_cast<long long>(closeness.size())) { Rcpp::stop("[ERROR]: index is not valid.\n"); } //

        //reference seqs all have indexes less than refEnd
    for (const long long it : closeness[index]) {
        if (!isRef[it]) closeSeqs.insert(it);  //you are a fit seq
    }

    return closeSeqs;
}
/***********************************************************************/
std::set<long long> OptiRefMatrix::getCloseRefSeqs(const long long index) {
    std::set<long long> closeSeqs;

    if (index < 0) return closeSeqs;
    if (index > static_cast<long long>(closeness.size())) { Rcpp::stop("[ERROR]: index is not valid.\n"); }
    //reference seqs all have indexes less than refEnd
    for (long long it : closeness[index]) {
        if (isRef[it]) { closeSeqs.insert(it); } //you are a ref seq
    }
    return closeSeqs;
}

std::map<std::string, long long> OptiRefMatrix::getNameIndexMap() const {

    const auto closenessSize = static_cast<long long>(closeness.size());
    std::map<std::string, long long> nameIndexes;
    for (long long i = 0; i < static_cast<long long>(nameMap.size()); i++) {
        std::vector<std::string> thisBinsSeqs;
        Utils::splitAtComma(nameMap[i], thisBinsSeqs);
        if (i < closenessSize) {
            nameIndexes[thisBinsSeqs[0]] = i;
        }
    }
    return nameIndexes;
}

std::unordered_set<long long> OptiRefMatrix::getIndexes(std::unordered_set<std::string> seqs) const {
    std::map<std::string, long long> nameIndexes = getNameIndexMap();
    std::map<std::string, long long>::iterator it;

    std::unordered_set<long long> indexes;
    for (const auto & seq : seqs) {
        it = nameIndexes.find(seq);
        if (it != nameIndexes.end()) {
            indexes.insert(it->second);
        }
    }
    return indexes;
}

/***********************************************************************/
//only used in open reference clustering
ListVector* OptiRefMatrix::getFitListSingle() {
    ListVector* singlelist = nullptr;

    if (singletons.empty()) return singlelist;
        singlelist = new ListVector();

    for (size_t i = 0; i < isSingleRef.size(); i++) {
        if (!isSingleRef[i]) {
            singlelist->push_back(singletons[i]); }
    }

    return singlelist;
}
/***********************************************************************/
void OptiRefMatrix::randomizeRefs() {
    const long long totalSeqs = (isRef.size()+isSingleRef.size());
    const long long numToSelect = totalSeqs * fitPercent;
    const auto refSingletonCutoff = static_cast<long long>(isRef.size());
    const auto singleSize = static_cast<long long>(isSingleRef.size());

    //select sequences to be reference
    std::set<long long> fitSeqsIndexes;
    if (!weights.empty()) {  fitSeqsIndexes = SubSample::getWeightedSample(weights, rng, numToSelect);  } //you have weighted selection
    else {
        long long numSelected = 0;
        while (numSelected < numToSelect) {
            fitSeqsIndexes.insert(Utils::getRandomIndex(totalSeqs - 1, rng)); //no repeats
            numSelected = static_cast<long long>(fitSeqsIndexes.size());
        }
    }

    //initilize isRef to true
    isRef.clear(); isRef.resize(refSingletonCutoff, true);
    isSingleRef.clear(); isSingleRef.resize(singleSize, true);

    //std::set isRef values
    for (const long long thisSeq : fitSeqsIndexes) {
        if (thisSeq < refSingletonCutoff) { //you are a non singleton seq in the closeness
            isRef[thisSeq] = false;
        }else { //thisSeq is a singleton
            isSingleRef[thisSeq-refSingletonCutoff] = false;
        }
    }

    //find number of fitDists, refDists and between dists
    calcCounts();
}
/***********************************************************************/
//for denovo method
int OptiRefMatrix::ReadFiles(const OptiData* matrix,
    const CountTableAdapter& adapter, std::unordered_set<std::string>& optionalRefNames) {
    std::set<long long> fitSeqsIndexes;
    long long count = 0;
    const std::vector<std::string> nameList = matrix->GetNameList();

    const auto numberOfSequences = static_cast<long long>(nameList.size());
    std::vector<double> abundances(numberOfSequences, 1);

    for (auto& name : nameList) {
        if (name.empty()) continue;
        if (refWeightMethod == "abundance")          { weights[count] = static_cast<long long>(adapter.FindTotalAbundance(name)); }
        else if (refWeightMethod == "connectivity")  { weights[count] = 1;          } //initialize
        else if (refWeightMethod == "accnos") { //fill fit indexes
            if (optionalRefNames.find(name) == optionalRefNames.end()) { //you are not a reference sequence
                fitSeqsIndexes.insert(count); //add as fit seq
            }
        }
        count++;
        nameMap.push_back(name);
    }


    singletons = matrix->GetSingletons();
    std::vector<bool> isSingletonVector(nameList.size(), true);
    const std::unordered_set<std::string> singletonNames = {singletons.cbegin(), singletons.cend()};
    for (size_t i = 0; i < isSingletonVector.size(); i++) {
        if (singletonNames.find(nameList[i]) != singletonNames.end()) continue;
        //if you are not a singleton
        isSingletonVector[i] = false;
    }

    closeness = matrix->GetCloseness();

    const long long numToSelect = static_cast<long long>(static_cast<float>(matrix->GetNameList().size()) * fitPercent);
    if (!weights.empty()) {
        fitSeqsIndexes = SubSample::getWeightedSample(weights, rng, numToSelect);
    } //you have weighted selection
    else {
        if (refWeightMethod != "accnos") { //fitIndexes are filled above { //randomly select references
            std::vector<long long> tempFitSeqs(numberOfSequences, 0);
            std::iota(tempFitSeqs.begin(), tempFitSeqs.end(), 0);
            Utils::mothurRandomShuffle(tempFitSeqs, rng);
            fitSeqsIndexes =  {tempFitSeqs.begin(), tempFitSeqs.begin() + numToSelect};
        }
    }


    //flag reference seqs singleton or not
    //Every empty in namemape is a singleton...Soo...

    for (long long i = 0; i < static_cast<long long>(nameList.size()); i++) {
        if (!isSingletonVector[i]) { // if you are not a singleton
            if (fitSeqsIndexes.find(i) != fitSeqsIndexes.end()) { //you are a fit seq
                isRef.push_back(false);
            }else {
                isRef.push_back(true);
            }
            continue;
        }
        // You are a singleton
        if (fitSeqsIndexes.find(i) != fitSeqsIndexes.end()) { //you are a fit seq singleton
            isSingleRef.push_back(false);
            continue;
        }
        isSingleRef.push_back(true);
    }

    calcCounts();

    return 0;
}
/***********************************************************************/
//for reading reference and fit files separately, reference method
int OptiRefMatrix::ReadFiles(const OptiData* refMatrix, const CountTableAdapter& refAdapter, const OptiData* fitMatrix,
    const CountTableAdapter& fitAdapter, const FastaDatabase& refFastaSequences,
    const FastaDatabase& fitFastaSequences) {
    std::unordered_map<std::string, long long> nameAssignment;

    std::vector<std::string> nameList = refAdapter.GetSequences();
    long long count = 0;
    for (const auto& name : nameList) {
        nameMap.emplace_back(name);
        nameAssignment[name] = count++;
    }

    long long refCount = count;
    singletons = refMatrix->GetSingletons();
    std::vector<bool> isSingletonVector(nameList.size(), true);
    std::unordered_set<std::string> singletonNames = {singletons.cbegin(), singletons.cend()};

    for (size_t i = 0; i < isSingletonVector.size(); i++) {
        if (singletonNames.find(nameList[i]) != singletonNames.end()) {
            continue;
        }
        //if you are not a singleton
        isSingletonVector[i] = false;
    }

    std::vector<std::string> fitNameList = fitAdapter.GetSequences();
    for (const auto& name : fitNameList) {
        // it->second = count; count++;
        nameMap.emplace_back(name);
        nameAssignment[name] = count++;
    }

    std::vector<std::string> fitSingletonNames = fitMatrix->GetSingletons();
    isSingletonVector.resize(count, true);
    singletonNames.insert(fitSingletonNames.begin(), fitSingletonNames.end());

    for (size_t i = 0; i < isSingletonVector.size(); i++) {
        if (singletonNames.find(nameList[i]) != singletonNames.end()) {
            continue;
        }
        //if you are not a singleton
        isSingletonVector[i] = false;
    }

    fitPercent = ((count-refCount) / static_cast<float>(count));

    //read bewtween file to update singletons
    // readColumnSingletons(singleton, betweendistfile, nameAssignment);
    closeness = refMatrix->GetCloseness();
    UpdateSingletons(refFastaSequences, fitFastaSequences, isSingletonVector, nameAssignment);

    long long nonSingletonCount = 0;
    std::map<long long, long long> singletonIndexSwap;
    for (long long i = 0; i < refCount; i++) {
        if (!isSingletonVector[i]) { //if you are not a singleton
            singletonIndexSwap[i] = nonSingletonCount;
            isRef.push_back(true);
            nonSingletonCount++;
        }else {
            singletons.push_back(nameMap[i]);
            isSingleRef.push_back(true);
        }
    }
    // refSingletonIndexSwap.clear();

    for (long long i = refCount; i < static_cast<long long>(isSingletonVector.size()); i++) {
        if (!isSingletonVector[i]) { //if you are not a singleton
            singletonIndexSwap[i] = nonSingletonCount;
            isRef.push_back(false);
            nonSingletonCount++;
        }else {
            singletons.push_back(nameMap[i]);
            isSingleRef.push_back(false);
        }
    }

    numSingletons = static_cast<long long>(singletons.size());
    std::map<std::string, std::string> names;
    calcCounts();

    return 0;
}

void OptiRefMatrix::calcCounts(){
    //find number of fitDists, refDists and between dists
    numRefDists = 0;
    numFitDists = 0;
    numBetweenDists = 0;
    numFitSingletons = 0;
    numFitSeqs = 0;
    numRefSingletons = 0;

    for (long long i = 0; i < static_cast<long long>(closeness.size()); i++) {

        const bool thisSeqIsRef = isRef[i];
        long long thisSeqsNumRefDists = 0;
        long long thisSeqsNumFitDists = 0;

        for (const long long& newB : closeness[i]) {
            if ((thisSeqIsRef) && (isRef[newB])) {  thisSeqsNumRefDists++; } //both refs
            else if ((thisSeqIsRef) && (!isRef[newB])) { numBetweenDists++; } // ref to fit dist
            else if ((!thisSeqIsRef) && (isRef[newB])) { numBetweenDists++; } // fit to ref dist
            else if ((!thisSeqIsRef) && (!isRef[newB])) { thisSeqsNumFitDists++; } // both fit
        }

        //a refSingleton or Fitsingleton may not be a true singleton (no valid dists in matrix), but may be a refSeq with no distances to other refs but distances to fitseqs. a fitsingleton may have dists to refs but no dists to other fitseqs.

        //you are a ref with no refdists, so you are a refsingleton
        if ((thisSeqIsRef) && (thisSeqsNumRefDists == 0)) {  numRefSingletons++; }
        else if ((!thisSeqIsRef) && (thisSeqsNumFitDists == 0)) {  numFitSingletons++; }
        else if ((!thisSeqIsRef) && (thisSeqsNumFitDists != 0)) {  numFitSeqs++; }

        numRefDists += thisSeqsNumRefDists;
        numFitDists += thisSeqsNumFitDists;
    }

    //counted twice
    numRefDists /= 2;
    numFitDists /= 2;
    numBetweenDists /= 2;
}
/***********************************************************************/


std::vector<long long> OptiRefMatrix::getCloseSeqs(const long long i){
    if (i < 0) {
        return {};
    }
    if (i > static_cast<long long>(closeness.size())) {
        return {};
    }
    return closeness[i];
}

std::string OptiRefMatrix::getName(const long long index) const {
    if (index < 0) {
        return "";
    }
    if (index > static_cast<long long>(closeness.size())) {
        return "";
    }
    return nameMap[index];
}

std::set<std::string> OptiRefMatrix::getNames(const std::unordered_set<long long> &indexes) const {
    std::set<std::string> names;
    for (const long long index : indexes) {
        names.insert(getName(index));
    }
    return names;
}
