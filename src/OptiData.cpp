//
// Created by Gregory Johnson on 3/29/24.
//

#include "MothurDependencies/OptiData.h"
bool OptiData::isClose(long long i, long long toFind){
    if (i < 0) { return false; }
    if (i > closeness.size()) { return false; }
    bool found = false;
    if (closeness[i].count(toFind) != 0) { found = true; }
    return found;

}
std::set<long long> OptiData::getCloseSeqs(long long i){
        if (i < 0) { std::set<long long> temp; return temp; }
        else if (i > closeness.size()) {  std::set<long long> temp; return temp; }
        return closeness[i];
}

long long OptiData::getNumClose(long long index)
{
    if (index < 0) { return 0; }
    else if (index > closeness.size()) { return 0; }
    else { return closeness[index].size(); }
}

ListVector* OptiData::getListSingle() {
    ListVector* singlelist = nullptr;
    if (singletons.size() == 0) { }
    else {
        singlelist = new ListVector();

        for (int i = 0; i < singletons.size(); i++) { singlelist->push_back(singletons[i]); }
    }
    return singlelist;

}

std::string OptiData::getName(long long index) {

    if (index < 0) { return ""; }
    else if (index > closeness.size()) { return ""; }

    return nameMap[index];

}