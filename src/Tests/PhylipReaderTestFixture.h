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
    bool TestReadPhylipFileFromRowData(const std::vector<RowData> &rowData, int expectedResult);

    bool TestReadPhylipFileToRowData(const CountTableAdapter& adapter,
        const std::string &filePath, const std::vector<RowData> &expectedResult);

    bool TestGetDistanceMatrix(const std::vector<RowData> &rowData, bool expectedResult);
    bool TestGetListVector(const std::vector<RowData> &rowData, int expectedResult);

private:
    void Setup() override;
    void TearDown() override;
    ReadPhylipMatrix* reader = nullptr;
};



#endif //PHYLIPREADERTESTFIXTURE_H
