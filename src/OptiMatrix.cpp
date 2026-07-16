
//
// Created by Gregory Johnson on 3/29/24.
//

#include "DataStructures/OptiMatrix.h"

bool OptiMatrix::isClose(const long long i, const long long toFind) const {
    if (i < 0) { return false; }
    if (i > static_cast<long long>(closeness.size())) {
        return false;
    }
    return std::binary_search(closeness[i].cbegin(), closeness[i].cend(), toFind);
   // return
    // if (closeness[i].find(toFind) != closeness[i].end()) { found = true; }
    // return found;

}
std::vector<long long> OptiMatrix::getCloseSeqs(const long long i) {
    if (i < 0) {
        std::vector<long long> temp; return temp;
    }
    if (i > static_cast<long long>(closeness.size())) {
        std::vector<long long> temp; return temp;
    }
        return closeness[i];
}

long long OptiMatrix::getNumClose(const long long index) const {
    if (index < 0)
        return 0;
    if (index > static_cast<long long>(closeness.size())) {
        return 0;
    }
    return closeness[index].size();

}

ListVector* OptiMatrix::getListSingle() const {
    ListVector* singlelist = nullptr;
    if (singletons.empty()) { }
    else {
        singlelist = new ListVector();

        for (const auto & singleton : singletons) {
            singlelist->push_back(singleton);
        }
    }
    return singlelist;

}

std::string OptiMatrix::getName(const size_t index) const {

    if (index > closeness.size()) {
        return "";
    }

    return nameMap.at(index);

}

bool OptiMatrix::mccValidCalc() const {

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
