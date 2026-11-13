//
// Created by Gregory Johnson on 11/3/26.
//

#ifndef REFACTOR_VECTORPAIRCOMPARE_H
#define REFACTOR_VECTORPAIRCOMPARE_H
#include <vector>

struct VectorPairCompare {
    bool operator()(const std::pair<long long, long long>& value, const long long& key) const {
        return (value.first < key);
    }
    bool operator()(const long long& key, const std::pair<long long, long long>& value) const {
        return (key < value.first);
    }

    bool operator()(const std::pair<std::string, long long>& value, const std::string& key) const {
        return (value.first < key);
    }
    bool operator()(const std::string& key, const std::pair<std::string, long long>& value) const {
        return (key < value.first);
    }
};
#endif //REFACTOR_VECTORPAIRCOMPARE_H