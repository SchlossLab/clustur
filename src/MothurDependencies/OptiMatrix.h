//
// Created by Gregory Johnson on 3/27/24.
//

#ifndef OPTIMATRIX_H
#define OPTIMATRIX_H
#include <vector>
#include<set>
#include <string>

#include "OptiData.h"


class OptiMatrix : public OptiData{

public:
    OptiMatrix(std::vector<std::set<long long>> close,  std::vector<std::string> name,
     std::vector<std::string> single, const double c): OptiData(c)
{
    closeness = close;
    nameMap = name;
    singletons = single;
}//closeness, namemap, singleton, cutoff
    OptiMatrix() = default;
    ~OptiMatrix(){}
    std::vector<std::set<long long>> GetCloseness() {return closeness;}
    std::vector<std::string> GetNameList() {return nameMap;}
    std::vector<std::string> GetSingletons() {return singletons;}
private:
    double cutoff{};
    std::string distFile, namefile, countfile, format, distFormat;
    bool sim{};
};



#endif //OPTIMATRIX_H
