//
// Created by Gregory Johnson on 8/30/24.
//

#include "Tests/DataFrameAdapterTestFixture.h"

bool DataFrameAdapterTestFixture::TestUnorderedMapToDataFrame(
    const std::unordered_map<std::string, std::vector<std::string>>& map,
    const std::vector<std::string>& expectedResult) {
    Setup();
    const Rcpp::DataFrame df = DataFrameAdapter::UnorderedMapToDataFrame(map);
    const std::vector<std::string> names = df.names();
    bool result = false;
    for(const auto& resultName : expectedResult) {
        result = std::find(names.begin(), names.end(), resultName) != names.end();
        if(!result)
            break;
    }
    TearDown();
    return result;
}

void DataFrameAdapterTestFixture::Setup() {
    // do nothing static method
}

void DataFrameAdapterTestFixture::TearDown() {
    // do nothing static method
}
