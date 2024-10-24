//
// Created by Gregory Johnson on 10/7/24.
//

#include "Adapters/DistanceFileReader.h"

Rcpp::DataFrame DistanceFileReader::SparseMatrixToDataFrame() const {
    const size_t size = rowDataMatrix.size();
    std::vector<std::string> indexOneNames(size * size);
    std::vector<std::string> indexTwoNames(size * size);
    std::vector<double> distances(size * size);

    size_t counter = 0;
    for(const auto& value : rowDataMatrix) {
        const std::string firstName = value.name;
        for(size_t i = 0; i < size; i++) {
            const double distance = value.rowValues[i];
            if (distance > cutoff || distance < 0) continue;
            indexOneNames[counter] = firstName;
            indexTwoNames[counter] = rowDataMatrix[i].name;
            distances[counter++] = distance;
        }
    }
    indexOneNames.resize(counter);
    indexTwoNames.resize(counter);
    distances.resize(counter);
    return Rcpp::DataFrame::create(Rcpp::Named("FirstName") = indexOneNames,
                                    Rcpp::Named("SecondName") = indexTwoNames,
                                    Rcpp::Named("Distance") = distances);
}

Rcpp::DataFrame DistanceFileReader::GetCountTable() const {
    // Recreate and return the count table
    // Might actually just store it in memory and return it
    // return countTable.GetCountTable();
    return countTable.GetCountTable();
}
void DistanceFileReader::SetRowDataMatrix(const std::vector<RowData>& data) {
    rowDataMatrix = data;
}
void DistanceFileReader::SetCountTable(CountTableAdapter data) {
    countTable = data;
}

bool DistanceFileReader::ReadRowDataMatrix(const std::vector<RowData> &rowData) {
    if (rowData.empty())
        return false;
    Utils util;
    std::string name = rowData[0].name;
    std::vector<std::string> matrixNames;
    const int nseqs = static_cast<int>(rowData.size());
    list = new ListVector(nseqs);
    list->set(0, name);

    sparseMatrix->resize(nseqs);
    for (int i = 1; i < nseqs; i++) {
        name = rowData[i].name;
        list->set(i, name);
        matrixNames.push_back(name); // Square matrix respresented as a sparse matrix. This causes a myraid of problems.
        // Where it should be giving me the distance between 1 and 2, its not giving me any distance, I removed them.
        for (int j = 0; j < i; j++) {
            auto distance = static_cast<float>(rowData[i].rowValues[j]);
            const bool equalivance = util.isEqual(distance, -1);
            if (equalivance) {
                distance = 1000000;
            } else if (sim) {
                distance = 1.0f - distance;
            }

            if (distance <= cutoff) {
                const PDistCell value(i, distance);
                sparseMatrix->addCell(j, value);
            }
        }
    }
    list->setLabel("0");
    return true;
}
