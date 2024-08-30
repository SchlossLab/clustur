//
// Created by Gregory Johnson on 8/30/24.
//

#ifndef DATAFRAMEADAPTERTESTFIXTURE_H
#define DATAFRAMEADAPTERTESTFIXTURE_H
#include <unordered_map>

#include "TestFixture.h"
#include "../Adapters/DataFrameAdapter.h"


class DataFrameAdapterTestFixture final : public TestFixture {
public:
    bool TestUnorderedMapToDataFrame(const std::unordered_map<std::string, std::vector<std::string>> &,
                                     const std::vector<std::string> &expectedResult);

private:
    void Setup() override;
    void TearDown() override;
    DataFrameAdapter* dataFrameAdapter = nullptr;
};



#endif //DATAFRAMEADAPTERTESTFIXTURE_H
