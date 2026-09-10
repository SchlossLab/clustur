//
// Created by Gregory Johnson on 9/10/26.
//

#ifndef REFACTOR_OPTIDATACOMPONENTFACTORY_H
#define REFACTOR_OPTIDATACOMPONENTFACTORY_H
#include "FastaDatabase.h"
#include "OptidataComponents.h"
#include "../MothurDependencies/PairwiseDistanceCalculator.h"

class OptiDataComponentFactory {
public:
    OptiDataComponentFactory() = default;
    ~OptiDataComponentFactory() = default;
    static OptiDataComponent CreateOptiDataComponent(const FastaDatabase& database, const PairwiseDistanceCalculator* calculator,
        const std::unordered_set<std::string>& filteredNameList) {
        const size_t filteredNameSize = filteredNameList.size();
        SparseDistanceMatrix* matrix = new SparseDistanceMatrix();
        ListVector* listVector = new ListVector();
        matrix->resize(filteredNameSize);
        listVector->resize(filteredNameSize);
        std::unordered_map<std::string, size_t> nameToIndex;
        nameToIndex.reserve(filteredNameSize);
        const std::vector<FastaData>& data = database.GetFastaDataBase();
        size_t count = 0;
        for (const auto& name : filteredNameList) {
            listVector->set(count, name);
            nameToIndex[name] = count++;
        }
        for (const auto& [name, sequence] : data) {
            if (filteredNameList.find(name) == filteredNameList.end()) continue;
            const size_t& index = nameToIndex[name];
            for (const auto& [otherName, otherSequence] : data) {
                if (filteredNameList.find(name) == filteredNameList.end()) continue;
                const size_t& otherIndex = nameToIndex[otherName];
                const double result = calculator->Execute(name, otherName);
                matrix->addCell(otherIndex, PDistCell(index, result));
            }
        }
        return {matrix, listVector};
    }
};
#endif //REFACTOR_OPTIDATACOMPONENTFACTORY_H