//
// Created by Gregory Johnson on 3/29/24.
//

#ifndef OPTIDATA_H
#define OPTIDATA_H
#include <vector>
#include <string>
#include <map>
#include <set>
#include <fstream>
#include <unordered_set>

#include "ListVector.h"


class OptiData {
public:
    OptiData(double c)  { cutoff = c; }
    OptiData() {}
    virtual ~OptiData(){}

    std::unordered_set<long long> getCloseSeqs(long long i);
    bool isClose(long long, long long) const;
    long long getNumClose(long long) const;
    std::string getName(long long); //name from nameMap index

    long long getNumSeqs() { return closeness.size(); }
    long long getNumSingletons() { return singletons.size(); }
    ListVector* getListSingle() const;
     virtual std::unordered_set<long long> getCloseRefSeqs(long long i) { std::unordered_set<long long> temp; return temp;  }
     virtual bool isCloseFit(long long j, long long i, bool&) { return false; }
protected:
    double cutoff;
    std::vector<std::unordered_set<long long>> closeness;  //closeness[0] contains indexes of seqs "close" to seq 0.
    std::vector<std::string> singletons; //name of seqs with NO distances in matrix, if name file is given then it contains 2nd column of namefile
    std::vector<std::string> nameMap;
};



#endif //OPTIDATA_H
