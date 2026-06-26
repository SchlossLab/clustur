//
// Created by Gregory Johnson on 6/15/26.
//

#ifndef REFACTOR_OPTIMATRIX_H
#define REFACTOR_OPTIMATRIX_H


#include <vector>
#include <string>
#include <unordered_set>
#include "ListVector.h"

class OptiMatrix final {

public:
    OptiMatrix(const std::vector<std::unordered_set<long long>>& close,  const std::vector<std::string>& name,
     const std::vector<std::string>& singleton, const double c)
    {
        closeness = close;
        nameMap = name;
        singletons = singleton;
        cutoff = c;
    }//closeness, namemap, singleton, cutoff
    std::vector<std::unordered_set<long long>> GetCloseness() {return closeness;}
    std::vector<std::string> GetNameList() {return nameMap;}
    std::vector<std::string> GetSingletons() {return singletons;}
    [[nodiscard]] size_t getNumSeqs() const { return closeness.size(); }
    [[nodiscard]] size_t getNumSingletons() const { return singletons.size(); }
    std::unordered_set<long long> getCloseSeqs(long long i);
    [[nodiscard]] bool isClose(long long, long long) const;
    [[nodiscard]] long long getNumClose(long long) const;
    [[nodiscard]] ListVector* getListSingle() const;
    [[nodiscard]] std::string getName(size_t) const; //name from nameMap index
    [[nodiscard]] bool isCloseFit(long long j, long long i, bool& k) const { return false; }
    bool mccValidCalc() const;
private:
    std::string distFile, namefile, countfile, format, distFormat;
    double cutoff;
    std::vector<std::unordered_set<long long>> closeness;  //closeness[0] contains indexes of seqs "close" to seq 0.
    std::vector<std::string> singletons; //name of seqs with NO distances in matrix, if name file is given then it contains 2nd column of namefile
    std::vector<std::string> nameMap;
};



#endif //REFACTOR_OPTIMATRIX_H