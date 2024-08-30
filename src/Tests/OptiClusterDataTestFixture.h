//
// Created by Gregory Johnson on 8/29/24.
//

#ifndef OPTICLUSTERDATATESTFIXTURE_H
#define OPTICLUSTERDATATESTFIXTURE_H
#include "TestFixture.h"
#include "../MothurDependencies/OpticlusterData.h"


class OptiClusterDataTestFixture final : public TestFixture {
public:
    bool CheckAddToData(OptiClusterInformation &information, bool expectedResult);
private:
    void Setup() override;
    void TearDown() override;
    OpticlusterData* data = nullptr;

};



#endif //OPTICLUSTERDATATESTFIXTURE_H
