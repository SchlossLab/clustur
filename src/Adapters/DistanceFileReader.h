//
// Created by Gregory Johnson on 10/7/24.
//

#ifndef DISTANCEFILEREADER_H
#define DISTANCEFILEREADER_H

#include <vector>

#include "CountTableAdapter.h"
#include "../MothurDependencies/ListVector.h"
#include "../RowData.h"
#include "../MothurDependencies/SparseDistanceMatrix.h"

class DistanceFileReader {
public:
    virtual ~DistanceFileReader() = default;
    // We need to deduce type, the easy way to do that is to see if there is a number for the first item read.
    // Phylip files have a number of sequences located at the top. We can use that to our advantage.
    virtual bool Read(const std::string& filePath) {return false;}
    DistanceFileReader(SparseDistanceMatrix*, ListVector*, bool);
    DistanceFileReader() = default;
    // Phylip files do not need a count table
    SparseDistanceMatrix* GetSparseMatrix() const {return new SparseDistanceMatrix(*sparseMatrix);}
    ListVector* GetListVector() const {return new ListVector(*list);}
    void SetRowDataMatrix(const std::vector<RowData>& data);
    void SetCountTable(CountTableAdapter data);
    Rcpp::DataFrame SparseMatrixToDataFrame() const;
    Rcpp::DataFrame GetCountTable() const;
    CountTableAdapter GetCountTableAdapter() const {return countTable;}
    double GetCutoff() const {return cutoff;}
    bool GetIsSimularity() const {return sim;}
    void CreateCountTableAdapter(const Rcpp::DataFrame&);



protected:
    SparseDistanceMatrix* sparseMatrix{};
    CountTableAdapter countTable;
    ListVector* list{};
    double cutoff = 0;
    bool sim = true;
private:


};

#endif //DISTANCEFILEREADER_H
