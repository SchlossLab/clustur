//
// Created by gregj on 8/27/2024.
//

#include "Adapters/CountTableAdapter.h"

void CountTableAdapter::CreateDataFrameMap(Rcpp::DataFrame countTable) {
    std::unordered_map<std::string, std::vector<double>> data;
    std::vector<std::string> columnNames = countTable.names();
    
    for(const auto& name : columnNames) {
        std::vector<double> columnData = countTable[name];
        data[name] = columnData;
    }
    dataFrameMap = data;
}
