//
// Created by Gregory Johnson on 9/16/24.
//

#ifndef COLUMNREADERTESTFIXTURE_H
#define COLUMNREADERTESTFIXTURE_H
#include "TestFixture.h"
#include <vector>
#include <string>
#include "../MothurDependencies/ColumnDistanceMatrixReader.h"


class ColumnReaderTestFixture final : TestFixture {
public:
    ColumnReaderTestFixture() = default;
    bool TestReadColumnFile(const std::string &file, const Rcpp::DataFrame& df, int expectedResult);
    bool TestGetDistanceMatrix(const std::string &filePath, const Rcpp::DataFrame& df, bool expectedResult);
    bool TestGetListVector(const std::string &filePath, const Rcpp::DataFrame& df, int expectedResult);

private:
    void Setup() override;
    void TearDown() override;
    ColumnDistanceMatrixReader* reader{};
};



#endif //COLUMNREADERTESTFIXTURE_H
