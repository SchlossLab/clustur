//
// Created by Gregory Johnson on 11/3/26.
//

#ifndef REFACTOR_VECTORPAIRCOMPARE_H
#define REFACTOR_VECTORPAIRCOMPARE_H
#include <vector>

struct VectorPairCompare {
    bool operator()(const std::pair<int, int>& value, const int& key) const {
        return (value.first < key);
    }
    bool operator()(const int& key, const std::pair<int, int>& value) const {
        return (key < value.first);
    }
};
#endif //REFACTOR_VECTORPAIRCOMPARE_H