//
// Created by Gregory Johnson on 8/29/24.
//

#include "Tests/SharedFileTestFixture.h"

bool SharedFileTestFixture::TestSharedFilePrintData(const std::vector<SharedAbundance> &data,
                                                    const Rcpp::DataFrame &expectedResult) {
    Setup();
    SharedFile sharedFile(data);
    Rcpp::DataFrame df = sharedFile.PrintData("otu");
    const std::vector<std::string> columnNames = df.names();
    const std::vector<std::string> expectedNames = expectedResult.names();
    TearDown();
    return columnNames == expectedNames;
}

void SharedFileTestFixture::Setup() {

}

void SharedFileTestFixture::TearDown() {

}
