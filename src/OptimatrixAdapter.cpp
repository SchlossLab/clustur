//
// Created by Gregory Johnson on 3/27/24.
//


// Recreate the countTable and Abund data files. We may need to do this in memory

#include "Adapters/OptimatrixAdapter.h"
#include "Adapters/CountTableAdapter.h"
#include "DataStructures/SparseDistanceMatrix.h"
#include "DataStructures/OptiMatrix.h"

OptiMatrix* OptimatrixAdapter::ConvertToOptimatrix(const SparseDistanceMatrix* matrixData,
                                                   const ListVector* listVector, const bool sim) const {
    const auto size = static_cast<long long>(matrixData->seqVec.size());
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
    std::vector<std::vector<long long>> closeness(nonSingletonCount);
    for(const auto& cell : matrixData->seqVec) {
        const std::string name = listVector->get(count + nameOffset);
        if(cell.empty()) {
            singletons.emplace_back(name);
            nameOffset++;
            continue;
        }
        std::vector<long long> cells(cell.size());
        size_t counter = 0;
        // cells.reserve(cell.size());
        for(const auto& row : cell) {
            float distance = row.dist;
            if (distance == -1) {
                distance = 1000000;
            } else if (sim) {
                distance = 1.0f - distance;
            }
            if(distance <= cutoff) {
                cells[counter++] = (indexSwap[row.index]);
                nameList[indexSwap[row.index]] = listVector->get(row.index);
            }
        }
        std::sort(cells.begin(), cells.end());
        closeness[count] = cells;
        count++;
    }
    return new OptiMatrix{closeness, nameList, singletons, cutoff};
}




