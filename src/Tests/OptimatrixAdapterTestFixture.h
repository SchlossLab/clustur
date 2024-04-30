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
    bool TestOptimatrixReturnsCorrectValue(const std::vector<int>&
        xPosition, const std::vector<int>& yPosition, const std::vector<double>& data);
    bool TestOptimatrixClosenessReturnsCorrectValue();
    bool TestOptimatrixSingletonReturnsCorrectValue();
    bool TestOptimatrixNameListReturnsCorrectValue();
private:
    void Setup() override;
    void TearDown() override;
    std::unique_ptr<OptimatrixAdapter> adapter = nullptr;
};



#endif //OPTIMATRIXTESTFIXTURE_H
