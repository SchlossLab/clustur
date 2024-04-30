//
// Created by gregj on 4/30/2024.
//

#ifndef OPTIMATRIXTESTFIXTURE_H
#define OPTIMATRIXTESTFIXTURE_H
#include "TestFixture.h"
#include "../Adapters/OptimatrixAdapter.h"

class OptimatrixAdapterTestFixture : public TestFixture{
public:
    void Setup() override;
    void TearDown() override;
private:
    OptimatrixAdapter* adapter = nullptr;
};



#endif //OPTIMATRIXTESTFIXTURE_H
