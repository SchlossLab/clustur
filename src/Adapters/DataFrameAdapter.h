//
// Created by Gregory Johnson on 8/29/24.
//

#ifndef DATAFRAMEADAPTER_H
#define DATAFRAMEADAPTER_H
#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <string>

class DataFrameAdapter {
public:
    static Rcpp::DataFrame UnorderedMapToDataFrame(const std::unordered_map<std::string, std::vector<std::string>> &map);
};



#endif //DATAFRAMEADAPTER_H
