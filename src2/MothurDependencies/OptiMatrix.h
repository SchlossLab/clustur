//
// Created by Gregory Johnson on 3/27/24.
//

#ifndef OPTIMATRIX_H
#define OPTIMATRIX_H
#include <vector>
#include<set>
#include <string>

#include "OptiData.h"


class OptiMatrix final : public OptiData{

public:
    OptiMatrix(const std::vector<std::unordered_set<long long>>& close,  const std::vector<std::string>& name,
     const std::vector<std::string>& singleton, const double c): OptiData(c)
{
    closeness = close;
    nameMap = name;
    singletons = singleton;
}//closeness, namemap, singleton, cutoff
    ~OptiMatrix() override = default;
    std::vector<std::unordered_set<long long>> GetCloseness() {return closeness;}
    std::vector<std::string> GetNameList() {return nameMap;}
    std::vector<std::string> GetSingletons() {return singletons;}
private:
    std::string distFile, namefile, countfile, format, distFormat;
};



#endif //OPTIMATRIX_H
