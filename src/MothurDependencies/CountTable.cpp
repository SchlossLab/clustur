//
// Created by gregj on 7/16/2024.
//

#include "CountTable.h"

int CountTable::readTable(const std::string &, const std::string &) {

    std::map<std::string, int> nameFileCounts = util.readNames(filename);
    for (std::map<std::string, int>::iterator it = nameFileCounts.begin(); it != nameFileCounts.end(); it++){
        indexNameMap[it->first] = uniques;
        totals.push_back(it->second);
        total += it->second;
        uniques++;
    }

}
