//
// Created by Gregory Johnson on 7/8/26.
//

#ifndef REFACTOR_OPTIREFMATRIX_H
#define REFACTOR_OPTIREFMATRIX_H

//
//  optirefmatrix.hpp
//  Mothur
//
//  Created by Sarah Westcott on 5/3/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//
#include <vector>
#include <string>
#include <set>
#include <unordered_set>
#include <map>

#include "FastaDatabase.h"
#include "ListVector.h"
#include "OptiData.h"


/* Looking to easily access ref, fit and combined information to compare OTU assignments for the references, the sequences to fit, and the merged reference fit OTUs */

class CountTableAdapter;

class OptiRefMatrix final : public OptiData {

public:
    OptiRefMatrix(const std::vector<std::vector<long long>>& close,  const std::vector<std::string>& name,
    const std::vector<std::string>& singleton, const double c);

    OptiRefMatrix(const OptiData* matrix, const CountTableAdapter& adapter, const double fP, const std::string& refWeight);
    OptiRefMatrix(const OptiData* matrix, const CountTableAdapter& adapter, std::unordered_set<std::string> accnosRefFileNames);

    OptiRefMatrix(const OptiData *refMatrix, const CountTableAdapter &refAdapter, const OptiData *fitMatrix,
                  const CountTableAdapter &fitAdapter, const FastaDatabase &refFastaSequences, const FastaDatabase &fitFastaSequences,
                  double cutoff);

    ~OptiRefMatrix() override = default;
    std::vector<long long> getTranslatedBins(std::vector<std::vector<std::string> >&, std::vector< std::vector<long long> >&) override;
    OptiData* extractMatrixSubset(std::unordered_set<long long>&);
    OptiData* extractMatrixSubset(std::unordered_set<std::string> &seqs);
    OptiData* extractRefMatrix();
    void randomizeRefs();
    std::vector<std::string> getRefSingletonNames();
    int ReadFiles(const OptiData *matrix,
                  const CountTableAdapter &adapter, std::unordered_set<std::string> &optionalRefNames);

    int ReadFiles(const OptiData *refMatrix, const CountTableAdapter &refAdapter, const OptiData *fitMatrix,
                  const CountTableAdapter &fitAdapter, const FastaDatabase &refFastaSequences, const FastaDatabase &fitFastaSequences);
    void UpdateSingletons(const FastaDatabase &refFastaSequences, const FastaDatabase &fitFastaSequences,
                          std::vector<bool> &singletons, std::unordered_map<std::string, long long> &nameMap);
    long long getNumFitTrueSingletons(); //reads that are true singletons (no valid dists in matrix) and are flagged as fit
    [[nodiscard]] long long getNumFitSingletons() const { return numFitSingletons; } //user singletons
    [[nodiscard]] long long getNumDists() const    { return (numFitDists +numRefDists+numBetweenDists); } //all distances under cutoff
    [[nodiscard]] long long getNumFitDists() const { return numFitDists; } //user distances under cutoff
    [[nodiscard]] long long getNumRefDists() const { return numRefDists; } //ref distances under cutoff
    std::unordered_set<long long> getIndexes(std::unordered_set<std::string> seqs) const;

    ListVector* getFitListSingle();

    std::vector<long long> getRefSeqs(); //every ref seq in matrix. Includes some that would have been singletons if not for the betweendistfile
    std::vector<long long> getFitSeqs(); //every fit seq in matrix. Includes some that would have been singletons if not for the betweendistfile
    [[nodiscard]] long long getNumFitSeqs() const { return numFitSeqs; } //only Fit seqs that are in fitdistfile and not singletons
    long long getNumFitClose(long long);
    long long getNumRefClose(long long);
    std::set<long long> getCloseFitSeqs(long long);
    std::set<long long> getCloseRefSeqs(long long);

    std::map<std::string, long long> getNameIndexMap() const;

    bool isCloseFit(long long, long long, bool&);
    std::vector<long long> getCloseSeqs(long long i);
    [[nodiscard]] std::string getName(long long) const; //name from nameMap index
    [[nodiscard]] std::set<std::string> getNames(const std::unordered_set<long long>& indexes) const;


protected:
    // double cutoff{};
    // std::vector<std::unordered_set<long long>> closeness;  //closeness[0] contains indexes of seqs "close" to seq 0.
    // std::vector<std::string> singletons; //name of seqs with NO distances in matrix, if name file is given then it contains 2nd column of namefile
    // std::vector<std::string> nameMap;
    std::map<long long, long long> weights; //seqeunce index in matrix to weight in chosing as reference
    std::string method, refWeightMethod;
    bool square;
    double cutoff;
    //a refSingleton or Fitsingleton may not be a true singleton (no valid dists in matrix), but may be a refSeq with no distances to other refs but distances to fitseqs. a fitsingleton may have dists to refs but no dists to other fitseqs.
    long long numFitDists, numRefDists, numRefSingletons, numFitSingletons, numBetweenDists, numSingletons, numFitSeqs;
    float fitPercent;
    std::vector<bool> isRef; //same size as closeness, this tells us whether a seq with distances in the matrix is a reference
    std::vector<bool> isSingleRef; ////same size as singletons, this tells us whether a seq WITHOUT distances in the matrix (singleton) is a reference
    void calcCounts();
};



#endif //REFACTOR_OPTIREFMATRIX_H