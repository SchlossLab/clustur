//
// Created by Gregory Johnson on 10/7/24.
//

#include <utility>

#include "FileReaders/DistanceFileReader.h"

#include "DataStructures/FastaDatabase.h"
#include "MothurDependencies/OneGapPairwiseDistance.h"
#include <RcppThread.h>
#include <mutex>


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

void DistanceFileReader::AddFittedDataToReference(
    const ListVector& otherListVector,
    const CountTableAdapter& otherCountTable,
    const FastaDatabase& database,
    const FastaDatabase& otherDatabase, const double cut, const int numberOfThreads) {
    const size_t otherSize = otherCountTable.GetSequences().size();
    const size_t currentSize = countTable.GetSequences().size();
    list.push_back(otherListVector);
    countTable.AddCountTable(otherCountTable);
    sparseMatrix.resize(currentSize + otherSize);
    const int newSize = otherSize + currentSize;
    list.resize(newSize);
    std::unordered_map<std::string, int> indexMap;
    indexMap.reserve(newSize);
    for (int i = 0; i < currentSize; i++) {
        indexMap[list.get(i)] = i;
    }
    PairwiseDistanceCalculator* calculator = new OneGapPairwiseDistance();
    int counter = currentSize;
    for (int i = 0; i < newSize; i++) {
        std::string names = otherListVector.get(i);
        std::vector<std::string> splitNames;
        Utils::splitAtComma(names, splitNames);
        for (const auto& splitName : splitNames) {
            list.set(counter, splitName);
            indexMap[splitName] = counter++;
        }
    }
    CalculateDistances(database, otherDatabase, calculator, indexMap, cut, numberOfThreads);
}
void DistanceFileReader::CalculateDistances(const FastaDatabase& database,
    const FastaDatabase& otherDatabase, const PairwiseDistanceCalculator* calculator,
    const std::unordered_map<std::string, int>& indexMap, const double cut, const int numberOfThreads) {
    const std::vector<FastaData>& refData = otherDatabase.GetFastaDataBase();
    const std::vector<FastaData>& fitData = database.GetFastaDataBase();
    std::mutex mutex;
    RcppThread::parallelFor(0, refData.size(), [&](size_t i) {
        const std::string& refName = refData[i].name;
        const std::string& refSequence = refData[i].sequence;
        const int iIndex = indexMap.at(refName);
        for (const auto &[fitName, fitSequence] : fitData) {
            if (const float result = static_cast<float>(calculator->Execute(refSequence,
                fitSequence)); result < cut) {
                const int jIndex = indexMap.at(fitName);
                mutex.lock();
                sparseMatrix.addCell(jIndex, {iIndex , result});
                mutex.unlock();
            }
        }
    }, numberOfThreads);
}
