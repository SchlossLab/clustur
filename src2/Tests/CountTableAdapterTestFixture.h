//
// Created by Gregory Johnson on 8/30/24.
//

#ifndef COUNTTABLEADAPTERTESTFIXTURE_H
#define COUNTTABLEADAPTERTESTFIXTURE_H
#include "TestFixture.h"
#include <string>
#include <Rcpp.h>
#include <vector>

#include "../Adapters/CountTableAdapter.h"

class CountTableAdapterTestFixture final : public TestFixture{
public:
    bool TestCreateDataFrameMap(const Rcpp::DataFrame &df, bool expectedResult);
    bool TestFindAbundanceBasedOnGroup(const std::string &groupName, const std::string &sampleName, const double expectedResult);
    bool TestFindAbundance(const std::string& sampleName, bool expectedResult);
    bool TestGetColumnByName(const std::string& columnName, bool expectedResult);
    bool TestGetGroups(const std::vector<std::string> &expectedResult);
    bool TestGetNameByIndex(int index, const std::string& expectedResult);
private:
    CountTableAdapter* countTable = nullptr;
    Rcpp::DataFrame dataframe;
    void Setup() override;
    void TearDown() override;
};



#endif //COUNTTABLEADAPTERTESTFIXTURE_H
