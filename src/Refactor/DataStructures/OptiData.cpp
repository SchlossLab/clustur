//
// Created by Gregory Johnson on 3/29/24.
//


#include "OptiData.h"

bool OptiData::isClose(const long long i, const long long toFind) const {
    if (i < 0) { return false; }
    if (i > static_cast<long long>(closeness.size())) {
        return false;
    }
    bool found = false;
    if (closeness[i].find(toFind) != closeness[i].end()) { found = true; }
    return found;

}
std::unordered_set<long long> OptiData::getCloseSeqs(const long long i){
    if (i < 0) {
        std::unordered_set<long long> temp; return temp;
    }
    if (i > static_cast<long long>(closeness.size())) {
        std::unordered_set<long long> temp; return temp;
    }
    return closeness[i];
}

size_t OptiData::getNumClose(const long long index) const {
    if (index < 0)
        return 0;
    if (index > static_cast<long long>(closeness.size())) {
        return 0;
    }
    return closeness[index].size();

}

ListVector OptiData::getListSingle() const {
    ListVector result;
    if (singletons.empty()) return result;

    for (const auto & singleton : singletons) {
        result.push_back(singleton);
    }
    return result;
}

std::string OptiData::getName(const long long index) {

    if (index < 0) {
        return "";
    }
    if (index > static_cast<long long>(closeness.size())) {
        return "";
    }

    return nameMap[index];

}

bool OptiData::mccValidCalc() const {

    const auto numSeqs = static_cast<double>(getNumSeqs());
    const double numDists = numSeqs * (numSeqs-1)/2;

    size_t totalClose = 0;
    //for each sequence (singletons removed on read)
    for (const auto & closenessData : closeness) {
        totalClose += closenessData.size();
    }

    // Inital setup of all singletons - badState <- TN == 0, FP == 0, FN == totalClose/2, TP = 0
    // Inital setup of one otu - badState <- TN == 0, FP == 0, FN == 0, TP = totalClose/2
    if ((numDists - (static_cast<double>(totalClose)/2)) == 0) {
        return false;
    }

    return true;
}
