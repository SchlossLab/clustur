//
// Created by Gregory Johnson on 3/27/24.
//


// Recreate the countTable and Abund data files. We may need to do this in memory
#include "Adapters/OptimatrixAdapter.h"

#include "Adapters/CountTableAdapter.h"
#include "MothurDependencies/SparseDistanceMatrix.h"

OptiMatrix* OptimatrixAdapter::ConvertToOptimatrix(const SparseDistanceMatrix* matrixData,
    const ListVector* listVector, const bool sim) {
    const auto size = static_cast<long long>(matrixData->seqVec.size());
    // std::unordered_map<long long, long long> singletonIndexSwap;
    std::vector<std::string> nameList(size);
    std::vector<std::string> singletons;
    std::unordered_map<unsigned long, int> indexSwap;
    int nonSingletonCount = 0;
    for(int i = 0; i < size; i++) {
        if(!matrixData->seqVec[i].empty()) {
            indexSwap[static_cast<unsigned long>(i)] = nonSingletonCount++;
        }
    }
    int count = 0;
    int nameOffset = 0;
    std::vector<std::unordered_set<long long>> closeness(nonSingletonCount);
    for(const auto& cell : matrixData->seqVec) {
        const std::string name = listVector->get(count + nameOffset);
        // nameList[count + nameOffset] = name;
        if(cell.empty()) {
            singletons.emplace_back(name);
            nameOffset++;
            continue;
        }
        std::unordered_set<long long> cells;
        for(const auto& row : cell) {
            float distance = row.dist;
            if (distance == -1) {
                distance = 1000000;
            } else if (sim) {
                distance = 1.0f - distance;
            }
            if(distance <= cutoff) {
                cells.insert(indexSwap[row.index]);
                nameList[row.index] = listVector->get(row.index);
            }
        }
        closeness[count] = cells;
        count++;
    }
    return new OptiMatrix{closeness, nameList, singletons, cutoff};
}




