//
// Created by Gregory Johnson on 10/7/24.
//

#include "Adapters/DistanceFileReader.h"

Rcpp::DataFrame DistanceFileReader::SparseMatrixToDataFrame() const {
    const size_t size = sparseMatrix.size();
    std::vector<std::string> indexOneNames(size * size);
    std::vector<std::string> indexTwoNames(size * size);
    std::vector<double> distances(size * size);

    size_t counter = 0;
    for(const auto& value : sparseMatrix) {
        const std::string firstName = value.name;
        for(size_t i = 0; i < size; i++) {
            indexOneNames[counter] = firstName;
            indexTwoNames[counter] = sparseMatrix[counter].name;
            distances[counter++] = value.rowValues[i];
        }
    }

    return Rcpp::DataFrame::create(Rcpp::Named("FirstName") = indexOneNames,
                                    Rcpp::Named("SecondName") = indexTwoNames,
                                    Rcpp::Named("Distance") = distances);
}

Rcpp::DataFrame DistanceFileReader::GetCountTable() const {
    // Recreate and return the count table
    // Might actually just store it in memory and return it
    // return countTable.GetCountTable();
    return {};
}
void DistanceFileReader::SetSparseMatrix(const std::vector<RowData>& data) {
    sparseMatrix = data;
}
void DistanceFileReader::SetCountTable(const CountTableAdapter& data) {
    countTable = data;
}
