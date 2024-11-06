//
// Created by Gregory Johnson on 3/29/24.
//

#include "MothurDependencies/OptiData.h"
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

long long OptiData::getNumClose(const long long index) const {
    if (index < 0)
        return 0;
    if (index > static_cast<long long>(closeness.size())) {
        return 0;
    }
    return closeness[index].size();

}

ListVector* OptiData::getListSingle() const {
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

std::string OptiData::getName(const long long index) {

    if (index < 0) {
        return "";
    }
    if (index > static_cast<long long>(closeness.size())) {
        return "";
    }

    return nameMap[index];

}
