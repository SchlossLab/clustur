//
// Created by gregj on 4/30/2024.
//

#ifndef OPTIMATRIXTESTFIXTURE_H
#define OPTIMATRIXTESTFIXTURE_H
#include <memory>

#include "TestFixture.h"
#include "../Adapters/OptimatrixAdapter.h"
#include <vector>

class OptimatrixAdapterTestFixture : public TestFixture{
public:
    bool TestOptimatrixReturnsCorrectValue(const std::vector<int>&xPosition,
        const std::vector<int>& yPosition, const std::vector<double>& data);
    bool TestOptimatrixClosenessReturnsCorrectValue(const std::vector<int> &xPosition,
        const std::vector<int> &yPosition, const std::vector<double> &data, int expectedSizeOfList);
    bool TestOptimatrixSingletonReturnsCorrectValue(const std::vector<int> &xPosition,
        const std::vector<int> &yPosition, const std::vector<double> &data, int expectedSizeOfList);
    bool TestOptimatrixNameListReturnsCorrectValue(const std::vector<int> &xPosition,
        const std::vector<int> &yPosition, const std::vector<double> &data, int expectedSizeOfList);

    bool ConvertToOptimatrixWithRowData(const std::vector<RowData> &matrixData, bool sim,
                                        const std::vector<std::set<long long>> &expectedResult);

private:
    void Setup() override;
    void TearDown() override;
    OptimatrixAdapter* adapter = nullptr;
};



#endif //OPTIMATRIXTESTFIXTURE_H
