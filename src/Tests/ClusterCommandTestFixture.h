//
// Created by Gregory Johnson on 5/1/24.
//

#ifndef CLUSTERCOMMANDTESTFIXTURE_H
#define CLUSTERCOMMANDTESTFIXTURE_H
#include <memory>
#include "TestFixture.h"
#include "../MothurDependencies/ClusterCommand.h"
#include "../MothurDependencies/OptiMatrix.h"



class ClusterCommandTestFixture : TestFixture{
public:
    bool TestRunOptiClusterReturnsCorrectData(OptiMatrix* mockMatrix,
        const std::string& expectedResult);
    bool TestSetIterationsWorksCorrectly(int iterations, bool expectResult);
private:
    ClusterCommand* clusterCommand = nullptr;
    void Setup() override;
    void TearDown() override;
};



#endif //CLUSTERCOMMANDTESTFIXTURE_H
