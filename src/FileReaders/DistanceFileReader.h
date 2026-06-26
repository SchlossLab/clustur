//
// Created by Gregory Johnson on 10/7/24.
//

#ifndef DISTANCEFILEREADER_H
#define DISTANCEFILEREADER_H

#include <vector>

#include "../Adapters/CountTableAdapter.h"
#include "../DataStructures/ListVector.h"
#include "../DataStructures/SparseDistanceMatrix.h"

class DistanceFileReader {
public:
    virtual ~DistanceFileReader() = default;
    DistanceFileReader() = default;
    // We need to deduce type, the easy way to do that is to see if there is a number for the first item read.
    // Phylip files have a number of sequences located at the top. We can use that to our advantage.
    virtual bool Read(const std::string& filePath) {return false;}
    virtual std::set<std::string> GetFailureParameters() {return failureParameters;};
    DistanceFileReader(const SparseDistanceMatrix&, const ListVector&, CountTableAdapter , double, bool);
    DistanceFileReader(const SparseDistanceMatrix&, const ListVector&, double, bool);
    explicit DistanceFileReader(CountTableAdapter);
    // Phylip files do not need a count table
    [[nodiscard]] SparseDistanceMatrix* GetSparseMatrix() const {return new SparseDistanceMatrix(sparseMatrix);}
    [[nodiscard]] ListVector* GetListVector() const {return new ListVector(list);}
    [[nodiscard]] Rcpp::DataFrame SparseMatrixToDataFrame() const;
    void SetCountTableAdapter(const CountTableAdapter& adapter);
    [[nodiscard]] Rcpp::DataFrame GetCountTable() const;
    [[nodiscard]] CountTableAdapter GetCountTableAdapter() const {return countTable;}
    [[nodiscard]] double GetCutoff() const {return cutoff;}
    [[nodiscard]] bool GetIsSimularity() const {return sim;}



protected:
    SparseDistanceMatrix sparseMatrix{};
    CountTableAdapter countTable;
    ListVector list{};
    double cutoff = 0;
    bool sim = true;
    std::set<std::string> failureParameters; //If we fail to read the file or anything else
private:


};

#endif //DISTANCEFILEREADER_H
