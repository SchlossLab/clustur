//
// Created by Gregory Johnson on 6/15/26.
//

#ifndef REFACTOR_OPTIMATRIX_H
#define REFACTOR_OPTIMATRIX_H


#include <vector>
#include <string>
#include <unordered_set>
#include "ListVector.h"
#include "OptiData.h"
#include "OptiRefMatrix.h"

class OptiMatrix final : public OptiData {

public:
    OptiMatrix(const std::vector<std::vector<long long>>& close,  const std::vector<std::string>& name,
     const std::vector<std::string>& singleton, const double c)
    {
        closeness = close;
        nameMap = name;
        singletons = singleton;
        cutoff = c;
    }//closeness, namemap, singleton, cutoff
    OptiMatrix() = default;
protected:
    std::string distFile, namefile, countfile, format, distFormat;
    // double cutoff;
    // std::vector<std::unordered_set<long long>> closeness;  //closeness[0] contains indexes of seqs "close" to seq 0.
    // std::vector<std::string> singletons; //name of seqs with NO distances in matrix, if name file is given then it contains 2nd column of namefile
    // std::vector<std::string> nameMap;
};



#endif //REFACTOR_OPTIMATRIX_H