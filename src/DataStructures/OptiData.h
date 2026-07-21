//
// Created by Gregory Johnson on 6/9/26.
//

#ifndef REFACTOR_OPTIDATA_H
#define REFACTOR_OPTIDATA_H
#include <vector>
#include <string>
#include <unordered_set>
#include "../DataStructures/ListVector.h"


class OptiData {
public:
    OptiData() = default;
    virtual ~OptiData() = default;
    void MoveData(OptiData* other) {
        cutoff = other->cutoff;
        closeness = std::move(other->closeness);
        singletons = std::move(other->singletons);
        nameMap = std::move(other->nameMap);
    }

protected:
    double cutoff{};
    std::vector<std::unordered_set<long long>> closeness;  //closeness[0] contains indexes of seqs "close" to seq 0.
    std::vector<std::string> singletons; //name of seqs with NO distances in matrix, if name file is given then it contains 2nd column of namefile
    std::vector<std::string> nameMap;
};




#endif //REFACTOR_OPTIDATA_H