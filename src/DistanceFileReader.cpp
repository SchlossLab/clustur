//
// Created by Gregory Johnson on 10/7/24.
//

#include "Adapters/DistanceFileReader.h"


DistanceFileReader::DistanceFileReader(SparseDistanceMatrix *sparseDistanceMatrix,
    ListVector *listVector):sparseMatrix(sparseDistanceMatrix), list(listVector) {}

Rcpp::DataFrame DistanceFileReader::SparseMatrixToDataFrame() const {
    const size_t size = sparseMatrix->seqVec.size();
    std::vector<std::string> indexOneNames;
    std::vector<std::string> indexTwoNames;
    // std::unordered_map<long long, std::unordered_set<long long>> duplicateRows;
    std::vector<double> distances;

    indexOneNames.reserve(size * size); // The max size it can be
    indexTwoNames.reserve(size * size);
    distances.reserve(size * size);
    long long count = 0;
    for(const auto& value : sparseMatrix->seqVec) {
        const std::string firstName = list->get(count);
        for(const auto& rowVal : value) {
            const auto rowIndex = static_cast<long long>(rowVal.index);
            // const auto duplicates = duplicateRows[count];
            // const auto biDirectionDuplicates = duplicateRows[rowIndex];
            // Remove the bidirectional values
            // if(duplicates.find(rowIndex) != duplicates.end() ||
            //     biDirectionDuplicates.find(count) != biDirectionDuplicates.end()) {
            //     continue;
            // }
            const double distance = rowVal.dist;
            const std::string secondName = list->get(rowIndex);
            indexOneNames.emplace_back(firstName);
            indexTwoNames.emplace_back(secondName);
            distances.emplace_back(distance);
            // duplicateRows[count].insert(rowIndex);
            // duplicateRows[rowIndex].insert(count);
        }
        count++;
    }
    return Rcpp::DataFrame::create(Rcpp::Named("FirstName") = indexOneNames,
                                    Rcpp::Named("SecondName") = indexTwoNames,
                                    Rcpp::Named("Distance") = distances);
}
void DistanceFileReader::CreateCountTableAdapter(const Rcpp::DataFrame& countTableDataFrame) {
    // countTable = CountTableAdapter();
    countTable.CreateDataFrameMap(countTableDataFrame);
}

Rcpp::DataFrame DistanceFileReader::GetCountTable() const {
    // Recreate and return the count table
    // Might actually just store it in memory and return it
    // return countTable.GetCountTable();
    return countTable.GetCountTable();
}