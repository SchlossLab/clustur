//
// Created by Gregory Johnson on 8/29/24.
//

#ifndef CLUSTERDATATESTFIXTURE_H
#define CLUSTERDATATESTFIXTURE_H
#include "TestFixture.h"
#include "../MothurDependencies/ClusterData.h"


class ClusterDataTestFixture final : public TestFixture {
public:
    bool CheckAddToData(ClusterInformation &information, bool expectedResult);

private:
    void Setup() override;
    void TearDown() override;
    ClusterData* data = nullptr;

};



#endif //CLUSTERDATATESTFIXTURE_H
