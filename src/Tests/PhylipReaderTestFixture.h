//
// Created by gregj on 7/9/2024.
//

#ifndef PHYLIPREADERTESTFIXTURE_H
#define PHYLIPREADERTESTFIXTURE_H
#include "TestFixture.h"
#include <string>
#include <vector>
#include "../RowData.h"
#include "../MothurDependencies/ReadPhylipMatrix.h"


class PhylipReaderTestFixture final: public TestFixture{
public:
    bool TestReadPhylipFile(const std::string &file, int expectedResult);
    bool TestGetDistanceMatrix(const std::string &file, bool expectedResult);
    bool TestGetListVector(const std::string &file, int expectedResult);

private:
    void Setup() override;
    void TearDown() override;
    ReadPhylipMatrix* reader = nullptr;
};



#endif //PHYLIPREADERTESTFIXTURE_H
