//
// Created by Gregory Johnson on 5/1/24.
//

#ifndef CLUSTERCOMMANDTESTFIXTURE_H
#define CLUSTERCOMMANDTESTFIXTURE_H
#include "TestFixture.h"
#include "../MothurDependencies/ClusterCommand.h"
#include "../MothurDependencies/OptiMatrix.h"


class ClusterCommandTestFixture : TestFixture{
public:
    bool TestRunOptiClusterReturnsCorrectData(const OptiMatrix* mockMatrix,
        const std::string& expectedResult);
    bool TestSetIterationsWorksCorrectly(int iterations, int expectedIterations);
private:
    std::unique_ptr<ClusterCommand> clusterCommand;
    void Setup() override;
    void TearDown() override;
};



#endif //CLUSTERCOMMANDTESTFIXTURE_H
