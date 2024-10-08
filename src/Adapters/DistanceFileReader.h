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
    virtual std::vector<RowData> ReadToRowData(const CountTableAdapter& adapter, const std::string& filePath) {return {};};
    // Phylip files do not need a count table
    virtual std::vector<RowData> ReadToRowData(const std::string& filePath) { return {};};
    virtual SparseDistanceMatrix* GetSparseMatrix() const = 0;
    virtual ListVector* GetListVector() const = 0;
    void SetRowDataMatrix(const std::vector<RowData>& data);
    void SetCountTable(CountTableAdapter data);

    Rcpp::DataFrame SparseMatrixToDataFrame() const;
    Rcpp::DataFrame GetCountTable() const;
    CountTableAdapter GetCountTableAdapter() const {return countTable;}
    double GetCutoff() const {return cutoff;}
    bool GetIsSimularity() const {return sim;}
    std::vector<RowData> GetRowDataMatrix() const {return rowDataMatrix;}

    bool ReadRowDataMatrix(const std::vector<RowData>& rowData);



protected:
    SparseDistanceMatrix* sparseMatrix{};
    ListVector* list{};
    double cutoff = 0;
    bool sim = true;
private:
    std::vector<RowData> rowDataMatrix;
    CountTableAdapter countTable;
};

#endif //DISTANCEFILEREADER_H
