//
// Created by Gregory Johnson on 9/10/26.
//

#ifndef REFACTOR_OPTIDATACOMPONENTFACTORY_H
#define REFACTOR_OPTIDATACOMPONENTFACTORY_H
#include <RcppThread.h>
#include <mutex>
#include <atomic>
#include "FastaDatabase.h"
#include "OptidataComponents.h"
#include "../MothurDependencies/PairwiseDistanceCalculator.h"

class OptiDataComponentFactory {
public:
    OptiDataComponentFactory() = default;
    ~OptiDataComponentFactory() = default;
    static OptiDataComponent CreateOptiDataComponent(const FastaDatabase& database, const PairwiseDistanceCalculator* calculator,
         const std::vector<std::string>& filteredNameList, const double cutoff, const int numberOfThreads) {
        const size_t filteredNameSize = filteredNameList.size();
        SparseDistanceMatrix matrix;
        ListVector listVector;
        matrix.resize(filteredNameSize);
        listVector.resize(filteredNameSize);
        std::unordered_map<std::string, size_t> nameToIndex;
        nameToIndex.reserve(filteredNameSize);
        const std::vector<FastaData>& data = database.GetFastaDataBase();
        size_t count = 0;
        for (const auto& name : filteredNameList) {

            listVector.set(count, name);
            nameToIndex[name] = count++;
        }

        const size_t dataSize = data.size();
        std::mutex mutex;
        RcppThread::parallelFor(0, dataSize, [&filteredNameList, &matrix, &data, &mutex, &calculator,
            &cutoff, &dataSize, &nameToIndex](size_t i) {

            const auto&[name, sequence] = data[i];
            if (!std::binary_search(filteredNameList.begin(), filteredNameList.end(), name)) return;
            const size_t& index = nameToIndex[name];

            for (size_t j = i + 1; j < dataSize; j++) {
                const auto&[otherName, otherSequence] = data[j];
                if (!std::binary_search(filteredNameList.begin(), filteredNameList.end(), otherName))
                    continue;

                const size_t& otherIndex = nameToIndex[otherName];
                const double result = calculator->Execute(sequence,
                  otherSequence);

                if (result > cutoff) continue;
                mutex.lock();
                matrix.addCell(otherIndex, PDistCell(index, result));
                mutex.unlock();
            }
        }, numberOfThreads);

        return {matrix, listVector};
    }
};
#endif //REFACTOR_OPTIDATACOMPONENTFACTORY_H