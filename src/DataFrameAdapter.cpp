//
// Created by Gregory Johnson on 8/29/24.
//

#include "Adapters/DataFrameAdapter.h"

Rcpp::DataFrame DataFrameAdapter::UnorderedMapToDataFrame(const std::unordered_map<std::string,
    std::vector<std::string>>& map) {
    Rcpp::DataFrame df;
    for(const auto& value : map) {
        df.push_back(value.second, value.first);
    }
    return df;
}
