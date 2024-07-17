//
// Created by gregj on 7/16/2024.
//

#ifndef COUNTTABLE_H
#define COUNTTABLE_H

#include <string>
#include <vector>
#include <map>

#include "ListVector.h"


class CountTable {
 public:

        CountTable() = default;
        ~CountTable() = default;

        //reads and creates smart enough to eliminate groups with zero counts
        int readTable(const std::string&, const std::string&); //filename, format - if format=fasta, read fasta file and create unique table. if format=name, use namefile to create count table
        ListVector getListVector();

    private:
        std::string filename;
        Utils util;
        bool hasGroups, isCompressed;
        int total, uniques;
        std::vector<std::string> groups;
        std::vector<int> totals;
        std::vector<int> totalGroups;
        std::map<std::string, int> indexNameMap; //maps seqName -> vector index in counts. seq1 -> 1 would mean seq1's counts are stored in counts[1].
        std::map<std::string, int> indexGroupMap;

        int find(int seq, int group, bool returnNext); //returns index of intPair for group passed in. If group is not present in seq, returns -1
        int getAbund(int seq, int group); //returns abundance of intPair for seq and group passed in. If group is not present in seq, returns 0
};



#endif //COUNTTABLE_H
