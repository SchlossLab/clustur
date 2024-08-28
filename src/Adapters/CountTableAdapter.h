//
// Created by gregj on 8/27/2024.
//

#ifndef COUNTTABLEADAPTER_H
#define COUNTTABLEADAPTER_H
#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_map>

class CountTableAdapter {
public:
    CountTableAdapter() = default;
    void CreateDataFrameMap(Rcpp::DataFrame countTable);
    std::unordered_map<std::string, std::vector<double>> GetDataFrameMap(){return dataFrameMap;}
private:
    std::unordered_map<std::string, std::vector<double>> dataFrameMap;
};



#endif //COUNTTABLEADAPTER_H
