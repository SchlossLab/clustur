//
// Created by Gregory Johnson on 7/8/26.
//

#include "DataStructures/OptiRefMatrix.h"
#include <Rcpp.h>

#include "Adapters/CountTableAdapter.h"
#include "DataStructures/OptiMatrix.h"
#include "RNG/SubSample.h"
//
//  optirefmatrix.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/3/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//


/***********************************************************************/
// OptiRefMatrix::OptiRefMatrix(string d, string nc, string f, string df, double c, string fit, string fitnc, string fitf, string fitdf, string betweend, string betweendf) : OptiData(c) {
//
//     string refdistfile, refnamefile, refcountfile, refformat, refdistformat, fitdistfile, fitnamefile, fitcountfile, fitformat, fitdistformat, betweendistfile, betweendistformat;
//
//     refdistfile = d; refdistformat = df; refformat = f; fitdistfile = fit; fitdistformat = fitdf; fitformat = fitf; betweendistfile = betweend; betweendistformat = betweendf;
//
//     numFitSingletons = 0;
//     numRefSingletons = 0;
//     numSingletons = 0;
//     numBetweenDists = 0;
//     numFitDists = 0;
//     numRefDists = 0;
//     numFitSeqs = 0;
//
//     fitPercent = 0;
//     refWeightMethod = "none";
//
//     square = false;
//
//     if (refformat == "name") { refnamefile = nc; refcountfile = ""; }
//     else if (refformat == "count") { refcountfile = nc; refnamefile = ""; }
//     else { refcountfile = ""; refnamefile = ""; }
//
//     if (fitformat == "name") { fitnamefile = fitnc; fitcountfile = ""; }
//     else if (fitformat == "count") { fitcountfile = fitnc; fitnamefile = ""; }
//     else { fitcountfile = ""; fitnamefile = ""; }
//
//     readFiles(refdistfile, refnamefile, refcountfile, refformat, refdistformat, fitdistfile, fitnamefile, fitcountfile, fitformat, fitdistformat, betweendistfile, betweendistformat);
// }
/***********************************************************************/
//Since we are extracting a subset of the seqs some reads that may not have been singletons
OptiData* OptiRefMatrix::extractRefMatrix() {
    std::set<long long> seqs; for (long long i = 0; i < isRef.size(); i++) { if (isRef[i]) { seqs.insert(i); } }

    std::vector<std::string> subsetNameMap;
    std::vector<std::string> subsetSingletons;
    std::vector< std::vector<long long> > subsetCloseness;
    std::map<long long, long long> thisNameMap;
    std::map<long long, long long> nonSingletonNameMap;
    std::vector<bool> singleton; singleton.resize(seqs.size(), true);
    int count = 0;

    for (std::set<long long>::iterator it = seqs.begin(); it != seqs.end(); it++) {
        long long seqNum = *it;
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
    for (long long i = 0; i < singleton.size(); i++) {
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

    for (int i = 0; i < isSingleRef.size(); i++) { if (isSingleRef[i]) { subsetSingletons.push_back(singletons[i]); } }

    auto* unfittedMatrix = new OptiRefMatrix(subsetCloseness, subsetNameMap, subsetSingletons, cutoff);

    return unfittedMatrix;
}
/***********************************************************************/
//given std::set of names, pull out their dists and create optimatrix
OptiData* OptiRefMatrix::extractMatrixSubset(std::unordered_set<std::string>& seqs) {
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
    for (long long i = 0; i < singleton.size(); i++) {
        if (!singleton[i]) { //if you are a singleton
            nonSingletonNameMap[i] = nonSingletonCount;
            nonSingletonCount++;
        }else { seqs.erase(nonSingletonNameMap[i]);  subsetSingletons.push_back(getName(nonSingletonNameMap[i])); } //remove from unfitted
    }
    singleton.clear();

    subsetCloseness.resize(nonSingletonCount);
    for (auto it = seqs.begin(); it != seqs.end(); it++) {


        long long seqNum = *it;

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

    return new OptiRefMatrix(subsetCloseness, subsetNameMap, subsetSingletons, cutoff);
}
/***********************************************************************/
std::vector<long long> OptiRefMatrix::getTranslatedBins(std::vector<std::vector<std::string> > & binNames, std::vector<std::vector<long long> > & fixedBins) {

    fixedBins.clear();

    std::map<std::string, long long> nameIndexes;
    std::set<std::string> unique;
    for (long long i = 0; i < nameMap.size(); i++) { //std::vector of std::string representing the sequences in the matrix from the name file.
        std::vector<std::string> thisSeqsReps; Utils::splitAtComma(nameMap[i], thisSeqsReps); //split redundant names
        if (i < closeness.size()) {  nameIndexes[thisSeqsReps[0]] = i;  } //this is a sequence with distances in the matrix
        if (thisSeqsReps.size() == 1) { //you are unique
            unique.insert(thisSeqsReps[0]);
        }
    }

    for (long long i = 0; i < singletons.size(); i++) {
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
bool OptiRefMatrix::isCloseFit(const long long i, const long long toFind, bool& isFit){
    if (i < 0) return false;
    if (i > closeness.size()) { Rcpp::stop("[ERROR]: index is not valid.\n");}

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
    std::vector<long long> refSeqsIndexes;
    for (long long i = 0; i < isRef.size(); i++) {
        if (isRef[i]) { refSeqsIndexes.push_back(i); }
    }
    return refSeqsIndexes;
}
/***********************************************************************/
std::vector<std::string> OptiRefMatrix::getRefSingletonNames() {
    std::vector<std::string> refSeqsNames;

    for (long long i = 0; i < isSingleRef.size(); i++) {
        if (isSingleRef[i]) { refSeqsNames.push_back(singletons[i]); }
    }

    return refSeqsNames;
}
/***********************************************************************/
std::vector<long long> OptiRefMatrix::getFitSeqs() {
    std::vector<long long> fitSeqsIndexes;
    for (long long i = 0; i < isRef.size(); i++) {
        if (!isRef[i]) { fitSeqsIndexes.push_back(i);  }
    }
    return fitSeqsIndexes;

}
/***********************************************************************/
long long OptiRefMatrix::getNumFitTrueSingletons() {
    return std::count(isSingleRef.begin(), isSingleRef.end(), false);
}
/***********************************************************************/
long long OptiRefMatrix::getNumFitClose(long long index) {
        long long numClose = 0;
        if (index < 0) return numClose;
        if (index > closeness.size()) Rcpp::stop("[ERROR]: index is not valid.\n");
            //reference seqs all have indexes less than refEnd
        for (const long long it : closeness[index]) {
            if (!isRef[it]) {  numClose++; } //you are a fit seq
        }


        return numClose;
}
/***********************************************************************/
long long OptiRefMatrix::getNumRefClose(long long index) {
        long long numClose = 0;

        if (index < 0) return numClose;
        if (index > closeness.size()) Rcpp::stop("[ERROR]: index is not valid.\n");
            //reference seqs all have indexes less than refEnd
        for (const long long it : closeness[index]) {
            if (isRef[it]) {  numClose++; } //you are a ref seq
        }

        return numClose;
}
/***********************************************************************/
std::set<long long> OptiRefMatrix::getCloseFitSeqs(long long index){
    std::set<long long> closeSeqs;

    if (index < 0) return closeSeqs;
    if (index > closeness.size()) { Rcpp::stop("[ERROR]: index is not valid.\n"); } //

        //reference seqs all have indexes less than refEnd
    for (const long long it : closeness[index]) {
        if (!isRef[it]) closeSeqs.insert(it);  //you are a fit seq
    }

    return closeSeqs;
}
/***********************************************************************/
std::set<long long> OptiRefMatrix::getCloseRefSeqs(long long index) {
    std::set<long long> closeSeqs;

    if (index < 0) return closeSeqs;
    if (index > closeness.size()) { Rcpp::stop("[ERROR]: index is not valid.\n"); }
    //reference seqs all have indexes less than refEnd
    for (long long it : closeness[index]) {
        if (isRef[it]) { closeSeqs.insert(it); } //you are a ref seq
    }
    return closeSeqs;
}

std::map<std::string, long long> OptiRefMatrix::getNameIndexMap() {
    std::map<std::string, long long> nameIndexes;
    for (int i = 0; i < nameMap.size(); i++) {
        std::vector<std::string> thisBinsSeqs; Utils::splitAtComma(nameMap[i], thisBinsSeqs);
        if (i < closeness.size()) { nameIndexes[thisBinsSeqs[0]] = i;  }
    }
    return nameIndexes;
}

std::unordered_set<long long> OptiRefMatrix::getIndexes(std::unordered_set<std::string> seqs) {
    std::map<std::string, long long> nameIndexes = getNameIndexMap();
    std::map<std::string, long long>::iterator it;

    std::set<long long> indexes;
    for (const auto & seq : seqs) {
        it = nameIndexes.find(seq);
        if (it != nameIndexes.end()) {
            indexes.insert(it->second);
        }
    }
}

/***********************************************************************/
//only used in open reference clustering
ListVector* OptiRefMatrix::getFitListSingle() {
    ListVector* singlelist = nullptr;

    if (singletons.empty()) return singlelist;
        singlelist = new ListVector();

    for (int i = 0; i < isSingleRef.size(); i++) {
        if (!isSingleRef[i]) {
            singlelist->push_back(singletons[i]); }
    }

    return singlelist;
}
/***********************************************************************/
void OptiRefMatrix::randomizeRefs() {
    const long long totalSeqs = (isRef.size()+isSingleRef.size());
    const long long numToSelect = totalSeqs * fitPercent;
    const long long refSingletonCutoff = isRef.size();
    const long long singleSize = isSingleRef.size();

    //select sequences to be reference
    std::set<long long> fitSeqsIndexes;
    if (!weights.empty()) {  fitSeqsIndexes = SubSample::getWeightedSample(weights, numToSelect);  } //you have weighted selection
    else {
        long long numSelected = 0;
        while (numSelected < numToSelect) {
            fitSeqsIndexes.insert(Utils::getRandomIndex(totalSeqs-1)); //no repeats
            numSelected = fitSeqsIndexes.size();
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
    // std::string namefile, countfile;
    // if (dupsFormat == "name") { namefile = dupsFile; countfile = ""; }
    // else if (dupsFormat == "count") { countfile = dupsFile; namefile = ""; }
    // else { countfile = ""; namefile = ""; }
    //
    // std::map<std::string, long long> nameAssignment;
    // if (namefile != "") { Utils::readNames(namefile, nameAssignment); }
    // else  {
    //     CountTable ct; ct.readTable(countfile, false, true);
    //     std::map<std::string, int> temp = ct.getNameMap();
    //     for (std::map<std::string, int>::iterator it = temp.begin(); it!= temp.end(); it++) {  nameAssignment[it->first] = it->second; }
    // }

    //select sequences to be reference

    // refWeightMethod = "abundance";
    closeness = matrix->GetCloseness();
    std::set<long long> fitSeqsIndexes;
    long long count = 0;
    const std::vector<std::string> nameList = adapter.GetSequences();
    // nameList.insert(nameList.end(), singletons.begin(), singletons.end());
    const long long numberOfSequences = nameList.size();
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
        // it->second = count; count++;
        nameMap.push_back(name);
        // nameAssignment[it->first] = it->second;
    }
    // for (std::map<std::string, long long>::iterator it = nameAssignment.begin(); it!= nameAssignment.end(); it++) {
    //     if (refWeightMethod == "abundance")          { weights[count] = it->second; }
    //     else if (refWeightMethod == "connectivity")  { weights[count] = 1;          } //initialize
    //     else if (refWeightMethod == "accnos") { //fill fit indexes
    //         if (optionalRefNames.count(it->first) == 0) { //you are not a reference sequence
    //             fitSeqsIndexes.insert(count); //add as fit seq
    //         }
    //     }
    //     it->second = count; count++;
    //     nameMap.push_back(it->first);
    //     nameAssignment[it->first] = it->second;
    // }

    //read file to find singletons
    // std::vector<bool> singleton;
    // singleton.resize(count, true);
    // std::map<long long, long long> singletonIndexSwap;
    //
    // if (distFormat == "column")        {  singletonIndexSwap = readColumnSingletons(singleton, distFile, nameAssignment);           }
    // else if (distFormat == "phylip")   {  singletonIndexSwap = readPhylipSingletons(singleton, distFile, count, nameAssignment);    }
    //


    singletons = matrix->GetSingletons();
    std::vector<bool> isSingletonVector(nameList.size(), true);
    const std::unordered_set<std::string> singletonNames = {singletons.cbegin(), singletons.cend()};
    int nonSingletonCount = 0;
    for (int i = 0; i < isSingletonVector.size(); i++) {
        if (singletonNames.find(nameList[i]) != singletonNames.end()) {
            // singletonIndexSwap[i] = nonSingletonCount;
            // nonSingletonCount++;
            continue;
        }
        //if you are not a singleton
        isSingletonVector[i] = false;
    }


    // std::map<std::string, std::string> names;
    // if (namefile != "") {
    //     //update names for reference
    //     util.readNames(namefile, names);
    //     for (int i = 0; i < numSingletons; i++) {
    //         singletons[i] = names[singletons[i]];
    //     }
    // }

    //read reference file distances
    // bool hasName = false;
    // if (namefile != "") { hasName = true; }
    // if (distFormat == "column")        {  readColumn(distFile, hasName, names, nameAssignment, singletonIndexSwap);     }
    // else if (distFormat == "phylip")   {  readPhylip(distFile, hasName, names, nameAssignment, singletonIndexSwap);     }


    //randomly select the "fit" seqs
    const long long numToSelect = static_cast<long long>(static_cast<float>(matrix->GetNameList().size()) * fitPercent);
    if (!weights.empty()) {  fitSeqsIndexes = SubSample::getWeightedSample(weights, numToSelect);  } //you have weighted selection
    else {
        if (refWeightMethod != "accnos") { //fitIndexes are filled above { //randomly select references
            long long numSelected = 0;
            // const long long totalSeqs = numberOfSequences;
            std::vector<long long> fitSeqsIndexes2(numberOfSequences, 0);
            std::iota(fitSeqsIndexes2.begin(), fitSeqsIndexes2.end(), 0);
            Utils::mothurRandomShuffle(fitSeqsIndexes2);
            fitSeqsIndexes =  {fitSeqsIndexes2.begin(), fitSeqsIndexes2.begin() + numToSelect};
            // while (numSelected < numToSelect) {
            //     fitSeqsIndexes.insert(Utils::getRandomIndex(numberOfSequences-1)); //no repeats
            //     numSelected = fitSeqsIndexes.size();
            // }
        }
    }

    //flag reference seqs singleton or not
    //Every empty in namemape is a singleton...Soo...

    for (long long i = 0; i < nameList.size(); i++) {
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
    // for (long long i = 0; i < static_cast<long long>(nameList2.size()); i++) {
    //     if (!nameList2[i].empty()) { //if you are not a singleton
    //
    //        //its a reference
    //     }else {
    //         if (fitSeqsIndexes.find(i) != fitSeqsIndexes.end()) { //you are a fit seq singleton
    //             isSingleRef.push_back(false);
    //         }else { isSingleRef.push_back(true); } //its a singleton reference
    //     }
    // }
    //find number of fitDists, refDists and between dists
    calcCounts();

    return 0;
}
/***********************************************************************/
//for reading reference and fit files separately, reference method
// int OptiRefMatrix::readFiles(std::string refdistfile, std::string refnamefile, std::string refcountfile, std::string refformat, std::string refdistformat, std::string fitdistfile, std::string fitnamefile, std::string fitcountfile, std::string fitformat, std::string fitdistformat, std::string betweendistfile, std::string betweendistformat){
//     try {
//         std::map<std::string, long long> nameAssignment;
//         if (refnamefile != "") { util.readNames(refnamefile, nameAssignment); }
//         else  {
//             CountTable ct; ct.readTable(refcountfile, false, true);
//             std::map<std::string, int> temp = ct.getNameMap();
//             for (std::map<std::string, int>::iterator it = temp.begin(); it!= temp.end(); it++) {  nameAssignment[it->first] = it->second; }
//         }
//
//         long long count = 0;
//         for (std::map<std::string, long long>::iterator it = nameAssignment.begin(); it!= nameAssignment.end(); it++) {
//             it->second = count; count++;
//             nameMap.push_back(it->first);
//             nameAssignment[it->first] = it->second;
//         }
//
//         long long refCount = count;
//         std::vector<bool> singleton; singleton.resize(count, true); //resize will only std::set new elements to true
//         std::map<long long, long long> refSingletonIndexSwap; //index into
//         if (refdistformat == "column")        {  refSingletonIndexSwap = readColumnSingletons(singleton, refdistfile, nameAssignment);          }
//         else if (refdistformat == "phylip")   {  refSingletonIndexSwap = readPhylipSingletons(singleton, refdistfile, count, nameAssignment);   }
//
//         //read fit file to find singletons
//         std::map<long long, long long> fitSingletonIndexSwap;
//         std::map<std::string, long long> fitnameAssignment;
//         if (fitnamefile != "") { util.readNames(fitnamefile, fitnameAssignment); }
//         else  {
//             CountTable ct; ct.readTable(fitcountfile, false, true);
//             std::map<std::string, int> temp = ct.getNameMap();
//             for (std::map<std::string, int>::iterator it = temp.begin(); it!= temp.end(); it++) {  fitnameAssignment[it->first] = it->second; }
//         }
//
//         for (std::map<std::string, long long>::iterator it = fitnameAssignment.begin(); it!= fitnameAssignment.end(); it++) {
//             it->second = count; count++;
//             nameMap.push_back(it->first);
//             nameAssignment[it->first] = it->second;
//         }
//
//         singleton.resize(count, true);
//         if (fitdistformat == "column")        {  fitSingletonIndexSwap = readColumnSingletons(singleton, fitdistfile, nameAssignment);          }
//         else if (fitdistformat == "phylip")   {  fitSingletonIndexSwap = readPhylipSingletons(singleton, fitdistfile, count, nameAssignment);   }
//
//         fitPercent = ((count-refCount) / (float) count);
//
//         //read bewtween file to update singletons
//         readColumnSingletons(singleton, betweendistfile, nameAssignment);
//
//         long long nonSingletonCount = 0;
//         std::map<long long, long long> singletonIndexSwap;
//         for (long long i = 0; i < refCount; i++) {
//             if (!singleton[i]) { //if you are not a singleton
//                 singletonIndexSwap[i] = nonSingletonCount;
//                 isRef.push_back(true);
//                 nonSingletonCount++;
//             }else {
//                 singletons.push_back(nameMap[i]);
//                 isSingleRef.push_back(true);
//             }
//         }
//         refSingletonIndexSwap.clear();
//
//         for (long long i = refCount; i < singleton.size(); i++) {
//             if (!singleton[i]) { //if you are not a singleton
//                 singletonIndexSwap[i] = nonSingletonCount;
//                 isRef.push_back(false);
//                 nonSingletonCount++;
//             }else {
//                 singletons.push_back(nameMap[i]);
//                 isSingleRef.push_back(false);
//             }
//         }
//         singleton.clear();
//         fitSingletonIndexSwap.clear();
//
//         numSingletons = singletons.size();
//         closeness.resize(nonSingletonCount);
//
//         std::map<std::string, std::string> names;
//         if (refnamefile != "") { util.readNames(refnamefile, names); }
//
//         if (fitnamefile != "") {
//             std::map<std::string, std::string> fitnames;
//             util.readNames(fitnamefile, fitnames);
//
//             names.insert(fitnames.begin(), fitnames.end()); //copy fit names into names
//         }
//
//         if ((fitnamefile != "") || (refnamefile != "")) {
//             for (int i = 0; i < singletons.size(); i++) {
//                 std::map<std::string, std::string>::iterator it = names.find(singletons[i]);
//                 if (it != names.end()) { //update singletons
//                     singletons[i] = it->second;
//                 }
//             }
//         }
//
//         //read reference file distances
//         bool refHasName = false;
//         if (refnamefile != "") { refHasName = true; }
//         if (refdistformat == "column")        {  readColumn(refdistfile, refHasName, names, nameAssignment, singletonIndexSwap);     }
//         else if (refdistformat == "phylip")   {  readPhylip(refdistfile, refHasName, names, nameAssignment, singletonIndexSwap);     }
//
//
//         //read fit distances
//         bool fitHasName = false;
//         if (fitnamefile != "") { fitHasName = true; }
//         if (fitdistformat == "column")        {  readColumn(fitdistfile, fitHasName, names, nameAssignment, singletonIndexSwap);     }
//         else if (fitdistformat == "phylip")   {  readPhylip(fitdistfile, fitHasName, names, nameAssignment, singletonIndexSwap);     }
//
//
//         //read in between distances
//         bool hasName = fitHasName;
//         if (!hasName && refHasName) { hasName = true; } //if either the ref or fit has a name file then std::set hasName
//         if (betweendistformat == "column")        {  readColumn(betweendistfile, hasName, names, nameAssignment, singletonIndexSwap);     }
//         else if (betweendistformat == "phylip")   {  readPhylip(betweendistfile, hasName, names, nameAssignment, singletonIndexSwap);     }
//
//         //find number of fitDists, refDists and between dists
//         calcCounts();
//
//         return 0;
//     }
//     catch(exception& e) {
//         m->errorOut(e, "OptiRefMatrix", "readFiles");
//         exit(1);
//     }
// }
// /***********************************************************************/
// std::map<long long, long long> OptiRefMatrix::readColumnSingletons(std::vector<bool>& singleton, std::string distFile, std::map<std::string, long long>& nameAssignment){
//     try {
//
//         ifstream fileHandle; util.openInputFile(distFile, fileHandle);
//
//         std::string firstName, secondName;
//         double distance;
//         std::map<long long, long long> singletonIndexSwap;
//
//         while(fileHandle){  //let's assume it's a triangular matrix...
//
//             fileHandle >> firstName; gobble(fileHandle);
//             fileHandle >> secondName; gobble(fileHandle);
//             fileHandle >> distance;	gobble(fileHandle); // get the row and column names and distance
//
//             if (m->getDebug()) { cout << firstName << '\t' << secondName << '\t' << distance << endl; }
//
//             if (m->getControl_pressed()) {  break; }
//
//             if (util.isEqual(distance,-1)) { distance = 1000000; }
//
//             if(distance <= cutoff){
//                 std::map<std::string,long long>::iterator itA = nameAssignment.find(firstName);
//                 std::map<std::string,long long>::iterator itB = nameAssignment.find(secondName);
//
//                 if(itA == nameAssignment.end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the name or count file, please correct\n"); exit(1);  }
//                 if(itB == nameAssignment.end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the name or count file, please correct\n"); exit(1);  }
//
//                 long long indexA = (itA->second);
//                 long long indexB = (itB->second);
//                 singleton[indexA] = false;
//                 singleton[indexB] = false;
//                 singletonIndexSwap[indexA] = indexA;
//                 singletonIndexSwap[indexB] = indexB;
//             }
//         }
//         fileHandle.close();
//
//         return singletonIndexSwap;
//     }
//     catch(exception& e) {
//         m->errorOut(e, "OptiRefMatrix", "readColumnSingletons");
//         exit(1);
//     }
// }
// /***********************************************************************/
//
// std::map<long long, long long> OptiRefMatrix::readPhylipSingletons(std::vector<bool>& singleton, std::string distFile, long long& count, std::map<std::string, long long>& nameAssignment){
//     try {
//         float distance;
//         long long nseqs;
//         std::string name;
//         std::map<long long, long long> singletonIndexSwap;
//
//         ifstream fileHandle;
//         std::string numTest;
//
//         util.openInputFile(distFile, fileHandle);
//         fileHandle >> numTest >> name;
//         nameMap.push_back(name);
//         singletonIndexSwap[0] = 0;
//         nameAssignment[name] = 0;
//
//         if (!util.isContainingOnlyDigits(numTest)) { m->mothurOut("[ERROR]: expected a number and got " + numTest + ", quitting.\n"); m->setControl_pressed(true); return singletonIndexSwap; }
//         else { convert(numTest, nseqs); }
//
//         //square test
//         char d;
//         while((d=fileHandle.get()) != EOF){
//             if(isalnum(d)){ square = true; fileHandle.putback(d); for(int i=0;i<nseqs;i++){ fileHandle >> distance;  } break; }
//             if(d == '\n'){ square = false; break; }
//         }
//
//         singleton.resize((count+nseqs), true);
//         if(square == 0){
//
//             for(long long i=1;i<nseqs;i++){
//                 if (m->getControl_pressed()) {  break; }
//
//                 fileHandle >> name; nameMap.push_back(name); singletonIndexSwap[i] = i;  nameAssignment[name] = i;
//
//                 for(long long j=0;j<i;j++){
//
//                     fileHandle >> distance;
//
//                     if (util.isEqual(distance,-1)) { distance = 1000000; }
//
//                     if(distance <= cutoff){
//                         singleton[i] = false;
//                         singleton[j] = false;
//                     }
//                 }
//             }
//         }else{
//             for(long long i=1;i<nseqs;i++){
//                 if (m->getControl_pressed()) {  break; }
//
//                 fileHandle >> name; nameMap.push_back(name); singletonIndexSwap[i] = i; nameAssignment[name] = i;
//
//                 for(long long j=0;j<nseqs;j++){
//                     fileHandle >> distance;
//
//                     if (util.isEqual(distance,-1)) { distance = 1000000; }
//
//                     if(distance <= cutoff && j < i){
//                         singleton[i] = false;
//                         singleton[j] = false;
//                     }
//                 }
//             }
//         }
//         fileHandle.close();
//
//         count += nseqs;
//
//         return singletonIndexSwap;
//     }
//     catch(exception& e) {
//         m->errorOut(e, "OptiRefMatrix", "readPhylipSingletons");
//         exit(1);
//     }
// }
// /***********************************************************************/
// int OptiRefMatrix::readPhylip(std::string distFile, bool hasName, std::map<std::string, std::string>& names, std::map<std::string, long long>& nameAssignment, std::map<long long, long long>& singletonIndexSwap){
//     try {
//         long long nseqs;
//         std::string name;
//         double distance;
//
//         ifstream in; std::string numTest;
//         util.openInputFile(distFile, in);
//
//         in >> numTest >> name;
//
//         if (hasName) { name = names[name]; } //redundant names
//         nameMap[singletonIndexSwap[0]] = name;
//
//
//         if (!util.isContainingOnlyDigits(numTest)) { m->mothurOut("[ERROR]: expected a number and got " + numTest + ", quitting.\n"); m->setControl_pressed(true); return 0; }
//         else { convert(numTest, nseqs); }
//
//         //square test
//         char d;
//         while((d=in.get()) != EOF){
//             if(isalnum(d)){ square = true; in.putback(d); for(int i=0;i<nseqs;i++){ in >> distance;  } break; }
//             if(d == '\n'){ square = false; break; }
//         }
//
//         std::string line = "";
//         if(!square){
//
//             for(long long i=1;i<nseqs;i++){
//
//                 if (m->getControl_pressed()) {  break; }
//
//                 in >> name; gobble(in);
//
//                 if (hasName) { name = names[name]; } //redundant names
//                 nameMap[singletonIndexSwap[i]] = name;
//
//                 for(long long j=0;j<i;j++){
//
//                     in >> distance; gobble(in);
//
//                     if (util.isEqual(distance,-1)) { distance = 1000000; }
//
//                     if(distance <= cutoff){
//                         if (refWeightMethod == "connectivity") { //count dists
//                             weights[i]++; weights[j]++;
//                         }
//                         long long newB = singletonIndexSwap[j];
//                         long long newA = singletonIndexSwap[i];
//                         closeness[newA].insert(newB);
//                         closeness[newB].insert(newA);
//                     }
//                 }
//             }
//         }else{
//             for(long long i=0;i<nseqs;i++){ in >> distance;  } gobble(in);
//
//             for(long long i=1;i<nseqs;i++){
//                 if (m->getControl_pressed()) {  break; }
//
//                 in >> name; gobble(in);
//
//                 if (hasName) { name = names[name]; } //redundant names
//                 nameMap[singletonIndexSwap[i]] = name;
//
//                 for(long long j=0;j<nseqs;j++){
//                     in >> distance; gobble(in);
//
//                     if (util.isEqual(distance,-1)) { distance = 1000000; }
//
//                     if(distance <= cutoff && j < i){
//                         if (refWeightMethod == "connectivity") { //count dists
//                             weights[i]++; weights[j]++;
//                         }
//                         long long newB = singletonIndexSwap[j];
//                         long long newA = singletonIndexSwap[i];
//                         closeness[newA].insert(newB);
//                         closeness[newB].insert(newA);
//                     }
//                 }
//             }
//         }
//         in.close();
//
//         return 0;
//     }
//     catch(exception& e) {
//         m->errorOut(e, "OptiRefMatrix", "readPhylip");
//         exit(1);
//     }
// }
// /***********************************************************************/
//
// int OptiRefMatrix::readColumn(std::string distFile, bool hasName, std::map<std::string, std::string>& names, std::map<std::string, long long>& nameAssignment, std::map<long long, long long>& singletonIndexSwap){
//     std::string firstName, secondName;
//     double distance;
//
//     ifstream in; util.openInputFile(distFile, in);
//
//     while(in){  //let's assume it's a triangular matrix...
//
//         in >> firstName; gobble(in);
//         in >> secondName; gobble(in);
//         in >> distance;	gobble(in); // get the row and column names and distance
//
//         if (m->getDebug()) { cout << firstName << '\t' << secondName << '\t' << distance << endl; }
//
//         if (m->getControl_pressed()) {  in.close();   return 0; }
//
//         if (util.isEqual(distance,-1)) { distance = 1000000; }
//
//         if(distance <= cutoff){
//             std::map<std::string,long long>::iterator itA = nameAssignment.find(firstName);
//             std::map<std::string,long long>::iterator itB = nameAssignment.find(secondName);
//
//             if(itA == nameAssignment.end()){  m->mothurOut("AAError: Sequence '" + firstName + "' was not found in the name or count file, please correct\n"); exit(1);  }
//             if(itB == nameAssignment.end()){  m->mothurOut("ABError: Sequence '" + secondName + "' was not found in the name or count file, please correct\n"); exit(1);  }
//
//             long long indexA = (itA->second);
//             long long indexB = (itB->second);
//
//             if (refWeightMethod == "connectivity") { //count dists
//                 weights[indexA]++; weights[indexB]++;
//             }
//
//             long long newB = singletonIndexSwap[indexB];
//             long long newA = singletonIndexSwap[indexA];
//             closeness[newA].insert(newB);
//             closeness[newB].insert(newA);
//
//             if (hasName) {
//                 std::map<std::string, std::string>::iterator itName1 = names.find(firstName);
//                 std::map<std::string, std::string>::iterator itName2 = names.find(secondName);
//
//                 if (itName1 != names.end()) { firstName = itName1->second;  } //redundant names
//                 if (itName2 != names.end()) { secondName = itName2->second;  } //redundant names
//             }
//
//             nameMap[newA] = firstName;
//             nameMap[newB] = secondName;
//         }
//     }
//     in.close();
//
//     return 1;
// }
/***********************************************************************/

void OptiRefMatrix::calcCounts(){
    //find number of fitDists, refDists and between dists
    numRefDists = 0;
    numFitDists = 0;
    numBetweenDists = 0;
    numFitSingletons = 0;
    numFitSeqs = 0;
    numRefSingletons = 0;

    for (long long i = 0; i < closeness.size(); i++) {

        const bool thisSeqIsRef = isRef[i];
        long long thisSeqsNumRefDists = 0;
        long long thisSeqsNumFitDists = 0;

        for (long long newB : closeness[i]) {
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

