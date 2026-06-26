//
// Created by Gregory Johnson on 6/9/26.
//

#ifndef REFACTOR_OPTIDATA_H
#define REFACTOR_OPTIDATA_H
#include <vector>
#include <string>
#include <unordered_set>
#include "../DataStructures/ListVector.h"


class OptiData final {
public:
    explicit OptiData(const double c)  { cutoff = c; }
    OptiData() = default;
    virtual ~OptiData() = default;

    std::unordered_set<long long> getCloseSeqs(long long i);
    [[nodiscard]] bool isClose(long long, long long) const;
    [[nodiscard]] size_t getNumClose(long long) const;
    std::string getName(long long); //name from nameMap index

    // in the case of all distances being below the cutoff, the TN and FP will always be 0
    // This is because nothing is considered "far apart". The mcc score will always be 0.
    // In this case we need to select a different calculator to cluster
    bool mccValidCalc() const;
    [[nodiscard]] size_t getNumSeqs() const { return closeness.size(); }
    [[nodiscard]] size_t getNumSingletons() const { return singletons.size(); }
    [[nodiscard]] ListVector getListSingle() const;
    std::unordered_set<long long> getCloseRefSeqs(long long i) { std::unordered_set<long long> temp; return temp;  }
    bool isCloseFit(long long j, long long i, bool&) { return false; }
protected:
    double cutoff{};
    std::vector<std::unordered_set<long long>> closeness;  //closeness[0] contains indexes of seqs "close" to seq 0.
    std::vector<std::string> singletons; //name of seqs with NO distances in matrix, if name file is given then it contains 2nd column of namefile
    std::vector<std::string> nameMap;
};




#endif //REFACTOR_OPTIDATA_H