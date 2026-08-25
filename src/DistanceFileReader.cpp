//
// Created by Gregory Johnson on 10/7/24.
//

#include <utility>

#include "FileReaders/DistanceFileReader.h"

#include "DataStructures/FastaDatabase.h"
#include "MothurDependencies/OneGapPairwiseDistance.h"


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

void DistanceFileReader::AddFittedDataToReference(const SparseDistanceMatrix* otherSparseMatrix,
    const ListVector* otherListVector,
    const CountTableAdapter& otherCountTable,
    const FastaDatabase& database,
    const FastaDatabase& otherDatabase, const double cut) {

    list.push_back(*otherListVector);
    countTable.AddCountTable(otherCountTable);
    sparseMatrix.addCells(*otherSparseMatrix);
    const int newSize = list.size();
    std::unordered_map<std::string, int> indexMap;
    indexMap.reserve(newSize);
    PairwiseDistanceCalculator* calculator = new OneGapPairwiseDistance();
    for (int i = 0; i < newSize; i++) {
        indexMap[list.get(i)] = i;
    }
    const std::vector<FastaData>& refData = database.GetFastaDataBase();
    const std::vector<FastaData>& otherData = otherDatabase.GetFastaDataBase();
    for (const auto &[refName, refSequence] : refData) {
        const int iIndex = indexMap[refName];
        for (const auto &[fitName, fitSequence] : otherData) {
            if (const float result = static_cast<float>(calculator->Execute(refSequence,
                fitSequence)); result < cut) {
                const int jIndex = indexMap[fitName];
                sparseMatrix.addCell(jIndex, {iIndex , result});
            }
        }
    }

}


void DistanceFileReader::AddInBetweenData(
   ListVector& otherListVector,
    CountTableAdapter& otherCountTable,
    const FastaDatabase& database,
    const FastaDatabase& otherDatabase, const double cut) {
    const size_t otherSize = otherCountTable.GetSequences().size();
    const size_t currentSize = countTable.GetSequences().size();
    SparseDistanceMatrix otherSparseMatrix;
    otherSparseMatrix.resize(otherSize);
    otherSparseMatrix.addCells(sparseMatrix);
    sparseMatrix = otherSparseMatrix;
    otherCountTable.AddCountTable(countTable);
    countTable = otherCountTable;

    const int newSize = otherSize + currentSize;

    otherListVector.push_back(list);
    list = ListVector(newSize);
    std::unordered_map<std::string, int> indexMap;
    indexMap.reserve(newSize);


    int counter = 0;
    for (int i = 0; i < newSize; i++) {
        std::string names = otherListVector.get(i);
        std::vector<std::string> splitNames;
        Utils::splitAtComma(names, splitNames);
        for (const auto& splitName : splitNames) {
            list.set(counter, splitName);
            indexMap[splitName] = counter++;
        }
    }

    PairwiseDistanceCalculator* calculator = new OneGapPairwiseDistance();
    const std::vector<FastaData>& refData = database.GetFastaDataBase();
    const std::vector<FastaData>& otherData = otherDatabase.GetFastaDataBase();
    for (const auto &[refName, refSequence] : refData) {
        const int iIndex = indexMap[refName];
        for (const auto &[fitName, fitSequence] : otherData) {
            if (const float result = static_cast<float>(calculator->Execute(refSequence,
                fitSequence)); result < cut) {
                const int jIndex = indexMap[fitName];
                sparseMatrix.addCell(jIndex, {iIndex , result});
                }
        }
    }

}
