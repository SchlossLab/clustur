//
// Created by Gregory Johnson on 8/30/24.
//

#include "Tests/CountTableAdapterTestFixture.h"

bool CountTableAdapterTestFixture::TestCreateDataFrameMap(const Rcpp::DataFrame& df, const bool expectedResult) {
    Setup();
    const bool result = countTable->CreateDataFrameMap(df);
    TearDown();
    return result == expectedResult;
}

bool CountTableAdapterTestFixture::TestFindAbundanceBasedOnGroup(const std::string& groupName,
    const std::string& sampleName,
    const double expectedResult) {
    Setup();
    countTable->CreateDataFrameMap(dataframe);
    const double abundance = countTable->FindAbundanceBasedOnGroup(groupName, sampleName);
    TearDown();
    return expectedResult == abundance;
}

bool CountTableAdapterTestFixture::TestFindAbundance(const std::string &sampleName, const bool expectedResult) {
    Setup();
    countTable->CreateDataFrameMap(dataframe);
    const double abundance = countTable->FindTotalAbundance(sampleName);
    TearDown();
    return expectedResult == (abundance != -1);
}

bool CountTableAdapterTestFixture::TestGetColumnByName(const std::string& columnName, const bool expectedResult) {
    Setup();
    countTable->CreateDataFrameMap(dataframe);
    const std::vector<double> columnData = countTable->GetColumnByName(columnName);
    TearDown();
    return expectedResult == !columnData.empty();

}

bool CountTableAdapterTestFixture::TestGetGroups(const std::vector<std::string>& expectedResult) {
    Setup();
    countTable->CreateDataFrameMap(dataframe);
    const auto groups = countTable->GetGroups();
    TearDown();
    return expectedResult == groups;
}

void CountTableAdapterTestFixture::Setup() {
    countTable = new CountTableAdapter();
    const std::vector<std::string> compounds{"1", "2", "3", "4", "5"};
    const std::vector<double> total{10, 20, 30, 40, 50};
    dataframe = Rcpp::DataFrame::create(
        Rcpp::Named("Representative Sequence") = compounds,
        Rcpp::Named("total") = total,
        Rcpp::Named("nogroup") = total);
}

void CountTableAdapterTestFixture::TearDown() {
    delete(countTable);
    dataframe = Rcpp::DataFrame();
}
