//
// Created by Gregory Johnson on 8/29/24.
//

#include "Adapters/DataFrameAdapter.h"

Rcpp::DataFrame DataFrameAdapter::UnorderedMapToDataFrame(const std::unordered_map<std::string,
    std::vector<std::string>>& map) {
    Rcpp::DataFrame df = Rcpp::DataFrame::create();
    for(const auto&[fst, snd] : map) {
        df.push_back(snd, fst);
    }
    return df;
}
