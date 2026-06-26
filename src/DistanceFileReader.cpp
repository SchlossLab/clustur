//
// Created by Gregory Johnson on 10/7/24.
//

#include <utility>

#include "FileReaders/DistanceFileReader.h"


DistanceFileReader::DistanceFileReader(const SparseDistanceMatrix& sparseDistanceMatrix,
                                       const ListVector& listVector, CountTableAdapter countTable, const double cutoff, const bool isSim):sparseMatrix(sparseDistanceMatrix),
countTable(std::move(countTable)), list(listVector), cutoff(cutoff), sim(isSim) {}

DistanceFileReader::DistanceFileReader(const SparseDistanceMatrix& sparseDistanceMatrix,
    const ListVector& listVector, const double cutoff, const bool isSim) :sparseMatrix(sparseDistanceMatrix),
list(listVector), cutoff(cutoff), sim(isSim) {}

DistanceFileReader::DistanceFileReader(CountTableAdapter countTableAdapter):
countTable(std::move(countTableAdapter)) {}

Rcpp::DataFrame DistanceFileReader::SparseMatrixToDataFrame() const {
    const size_t size = sparseMatrix.seqVec.size();
    std::vector<std::string> indexOneNames;
    std::vector<std::string> indexTwoNames;
    std::vector<double> distances;

    std::vector<bool> hasComputedRowDistances(size, false);
    indexOneNames.reserve(size * size); // The max size it can be
    indexTwoNames.reserve(size * size);
    distances.reserve(size * size);
    long long count = 0;
    for(const auto& value : sparseMatrix.seqVec) {
        const std::string firstName = list.get(count);
        for(const auto& rowVal : value) {
            const auto rowIndex = static_cast<long long>(rowVal.index);
            if(hasComputedRowDistances[rowIndex])
                continue;
            const double distance = rowVal.dist;
            const std::string secondName = list.get(rowIndex);
            indexOneNames.emplace_back(firstName);
            indexTwoNames.emplace_back(secondName);
            distances.emplace_back(distance);
        }
        hasComputedRowDistances[count++] = true;
    }
    return Rcpp::DataFrame::create(Rcpp::Named("FirstName") = indexOneNames,
                                    Rcpp::Named("SecondName") = indexTwoNames,
                                    Rcpp::Named("Distance") = distances);
}

void DistanceFileReader::SetCountTableAdapter(const CountTableAdapter &adapter) {
    countTable = adapter;
}

Rcpp::DataFrame DistanceFileReader::GetCountTable() const {
    return countTable.ReCreateDataFrame();
}
