//
// Created by Gregory Johnson on 3/27/24.
//


// Recreate the countTable and Abund data files. We may need to do this in memory
#include "Adapters/OptimatrixAdapter.h"

#include "Adapters/CountTableAdapter.h"
#include "MothurDependencies/CountTable.h"
#include "MothurDependencies/SparseDistanceMatrix.h"

//TODO: We have to create a full matrix, and add the values as from the sparse matrix inside the new matrix.
//TODO: This will create the full matrix, then and only then, will we be able to create a list of closeness,
//TODO: singletons, and names properly.
//Ensure we are not using a square matrix
//Fix the singletonIssue
//
// Hard code for matrix
// WE are going to add test cases in
//DONE Add TESTING SUITE FOR EASY DEBUGGING AND TESTING IN PURE CPP
OptiMatrix *OptimatrixAdapter::ConvertToOptimatrix(const std::vector<int> &xPosition,
                                                                    const std::vector<int> &yPosition,
                                                                    const std::vector<double> &data) {
    // The closeness map does not contain sequences that have 0 values....
    //If there are no values within the range of it, and its empty, REMOVE IT
    // All the values that are correlated to this also belong in the nameMap...
    // This was the issue!
    // Create the singleton vector
    // Create the closeness vector
    // Create the names vector
    // Creating these will allow us to have a prototype opticluster sparse matrix.
    const size_t count = data.size();
    // Get unique names
    // TYPE CONVERSIONS TO STRINGS ARE SLOW, but for a one-to-one concept, we are converting to strings
    // As per the OptiMatrix class
    std::map<int, std::set<long long>> closenessMap;
    std::vector<std::string> singletonList;
    // Shouldn't be duplicates, but there can be duplicates. There is a way to deal with it
    std::unordered_set<int> singletonCandidates; // If your names is not included in this set you are an singleton
    std::set<int> names;
    for (size_t i = 0; i < count; i++) {
        int currentXPos = xPosition[i];
        int currentYPos = yPosition[i];
        if (data[i] > cutoff) // its 1 - score, so the higher your score, the closer to the cutoff you are
        {
            // You are a singleton, so we will add you to the list
            //singletonList.emplace_back(std::to_string(xPosition[i]));
            // At this position, I dont cluster with the current x,y Position, so I am a singleton
            //List of possible singletons
            //I am using the index not the actually position, the xPos is the position to the other.
            singletonCandidates.insert(currentXPos);
            singletonCandidates.insert(currentYPos);
            continue;
        }
        // xPosition[i] is the name in this context
        // Linked list?
        closenessMap[currentXPos].emplace(currentYPos);
        closenessMap[currentYPos].emplace(currentXPos);
    }
    //TODO Change this back into a vector, we do not need to delete values
    std::vector<std::unordered_set<long long> > adjustedClosenessList;
    for (const auto &closenessValues: closenessMap) {
        adjustedClosenessList.emplace_back(closenessValues.second.begin(), closenessValues.second.end());
        names.insert(closenessValues.first);
    }
    //TODO Work on creating singletons properly!
    for (const auto &key: singletonCandidates) {
        if(names.find(key) != names.end()) //If its not found, its a singleton
            continue;
        singletonList.emplace_back(std::to_string(key));
    }
    int counter = 0;
    //Convert the unique singletons over to a list
    std::vector<std::string> namesVector(names.size());
    for (const auto &name: names) {
        namesVector[counter++] = std::to_string(name);
    }
    std::vector<std::unordered_set<long long> > closeness = adjustedClosenessList;
    std::vector<std::string> nameList = namesVector;
    std::vector<std::string> singletons = singletonList;
    return new OptiMatrix{adjustedClosenessList, namesVector, singletonList, cutoff};
}
// TODO Change RowData to sparse matix
OptiMatrix* OptimatrixAdapter::ConvertToOptimatrix(const std::vector<RowData>& matrixData, const bool sim) {
    const auto size = static_cast<long long>(matrixData.size());
    std::vector<std::unordered_set<long long>> closeness;
    std::vector<std::string> nameList;
    std::vector<std::string> singletons;
    Utils util;
    std::vector<bool> singletonList(size, true);
    nameList.resize(size);
    std::unordered_map<long long, long long> singletonIndexSwap;
    for(long long i = 0; i < size; i++) {
        nameList[i] = matrixData[i].name;
        singletonIndexSwap[i] = i;
        for(long long j = 0; j < i; j++) {
            auto distance = static_cast<float>(matrixData[i].rowValues[j]);
            const bool equalivance = util.isEqual(distance, -1);
            if (equalivance) {
                distance = 1000000;
            } else if (sim) {
                distance = 1.0f - distance;
            }
            if(distance <= cutoff) {
                singletonList[i] = false; // Find out who is a singleton
                singletonList[j] = false;
                singletonIndexSwap[i] = i;
                singletonIndexSwap[j] = j;
            }

        }
    }
    int nonSingletonCount = 0;
    for(size_t i = 0; i < singletonList.size(); i ++) {
        if(!singletonList[i]) {
            singletonIndexSwap[static_cast<long long>(i)] = nonSingletonCount;
            nonSingletonCount++;
        } //Remove all singletonss
        else
            singletons.emplace_back(matrixData[i].name);
    }
    closeness.resize(nonSingletonCount);
    for(long long i = 0; i < size; i++) {
        nameList[singletonIndexSwap[i]] = matrixData[i].name;
        for(long long j = 0; j < i; j++) {
            auto distance = static_cast<float>(matrixData[i].rowValues[j]);
            const bool equalivance = util.isEqual(distance, -1);
            if (equalivance) {
                distance = 1000000;
            } else if (sim) {
                distance = 1.0f - distance;
            }
            if(distance <= cutoff) {
                long long newB = singletonIndexSwap[j];
                long long newA = singletonIndexSwap[i];
                closeness[newA].insert(newB);
                closeness[newB].insert(newA);
            }

        }
    }
    return new OptiMatrix{closeness, nameList, singletons, cutoff};
}

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
    std::vector<std::unordered_set<long long>> closeness(nonSingletonCount);
    for(const auto& cell : matrixData->seqVec) {
         const std::string name =listVector->get(count);
         nameList[count] = name;
        if(cell.empty()) {
            singletons.emplace_back(name);
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
            }
        }
        closeness[count] = cells;
        count++;
    }



















    // Utils util;
    // std::vector<bool> singletonList(size, true);
    // nameList.resize(size);
    // std::unordered_map<long long, long long> singletonIndexSwap;
    // for(long long i = 0; i < size; i++) {
    //     nameList[i] = listVector->get(i);
    //     singletonIndexSwap[i] = i;
    //     for(long long j = 0; j < static_cast<long long>(matrixData->seqVec[i].size()); j++) {
    //         const auto cell = matrixData->seqVec[i][j];
    //         auto distance = static_cast<float>(cell.dist);
    //         const bool equalivance = util.isEqual(distance, -1);
    //         if (equalivance) {
    //             distance = 1000000;
    //         } else if (sim) {
    //             distance = 1.0f - distance;
    //         }
    //         if(distance <= cutoff) {
    //             singletonList[i] = false; // Find out who is a singleton
    //             singletonList[cell.index] = false;
    //             singletonIndexSwap[i] = i;
    //             singletonIndexSwap[static_cast<long long>(cell.index)] = static_cast<long long>(cell.index);
    //         }
    //
    //     }
    // }
    // int nonSingletonCount = 0;
    // for(size_t i = 0; i < singletonList.size(); i ++) {
    //     if(!singletonList[i]) {
    //         singletonIndexSwap[static_cast<long long>(i)] = nonSingletonCount++;
    //     } //Remove all singletonss
    //     else
    //         singletons.emplace_back(listVector->get(static_cast<long long>(i)));
    // }
    // closeness.resize(nonSingletonCount);
    // for(long long i = 0; i < size; i++) {
    //     nameList[singletonIndexSwap[i]] = listVector->get(i);
    //     for(long long j = 0; j < static_cast<long long>(matrixData->seqVec[i].size()); j++) {
    //         const PDistCell cell = matrixData->seqVec[i][j];
    //         auto distance = static_cast<float>(cell.dist);
    //         const bool equalivance = util.isEqual(distance, -1);
    //         if (equalivance) {
    //             distance = 1000000;
    //         } else if (sim) {
    //             distance = 1.0f - distance;
    //         }
    //         if(distance <= cutoff) {
    //             long long newB = singletonIndexSwap[static_cast<long long>(cell.index)];
    //             long long newA = singletonIndexSwap[i];
    //             closeness[newA].insert(newB);
    //             closeness[newB].insert(newA);
    //         }
    //
    //     }
    // }
    return new OptiMatrix{closeness, nameList, singletons, cutoff};
}




