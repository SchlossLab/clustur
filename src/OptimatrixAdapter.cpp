//
// Created by Gregory Johnson on 3/27/24.
//


// Recreate the countTable and Abund data files. We may need to do this in memory
#include "Adapters/OptimatrixAdapter.h"

//TODO
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
    std::vector<std::pair<int, std::set<long long> > > closenessList(count);

    std::vector<std::string> singletonList;
    // Shouldn't be duplicates, but there can be duplicates. There is a way to deal with it
    std::unordered_map<int, std::vector<int> > singletonClusterMap;
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
            singletonClusterMap[currentXPos].emplace_back(currentYPos);
            continue;
        }
        // xPosition[i] is the name in this context
        // Linked list?
        closenessList[currentXPos].first = currentXPos;
        closenessList[currentXPos].second.insert(yPosition[i]);
    }
    //TODO Change this back into a vector, we do not need to delete values
    const std::list<std::pair<int, std::set<long long>>> closenessListWithoutEmpties(closenessList.begin(),
        closenessList.end());
    std::vector<std::set<long long> > adjustedClosenessList;
    closenessList.clear();
    for (auto &[index, closenessValues]: closenessListWithoutEmpties) {
        if (!closenessValues.empty()) {
            adjustedClosenessList.emplace_back(closenessValues);
            names.insert(index);
            continue;
        }
    }

    const size_t maxSizeOfCluster = names.size();
    for (const auto &[key, value]: singletonClusterMap) {
        if (value.size() < maxSizeOfCluster) {
            // If the number of values that you are a singleton to is greater
            //than the amount of currently avaliable names (Means you are completely alone), than you are a singleton
            //Means I am not a singleton
            continue;
        }
        singletonList.emplace_back(std::to_string(key));
    }
    int counter = 0;
    //Convert the unique singletons over to a list
    std::vector<std::string> namesVector(names.size());
    for (const auto &name: names) {
        namesVector[counter++] = std::to_string(name);
    }
    closeness = adjustedClosenessList;
    nameList = namesVector;
    singletons = singletonList;
    return new OptiMatrix{adjustedClosenessList, namesVector, singletonList, cutoff};
}