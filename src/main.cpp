#include <iostream>
#include <fstream>
//TODO Next week: We have to separate the code out into a package! We have all the source code, now need to consider how
//it will look in a package!
#include "Adapters/OptimatrixAdapter.h"
#include "MothurDependencies/ListVector.h"
#include "TestHelpers/TestHelper.h"
#include "Tests/OptimatrixAdapterTestFixture.h"
#if DEBUG_RCPP
#include <vector>
#include "Adapters/OptimatrixAdapter.h"
#include "MothurDependencies/OptiMatrix.h"
#include "MothurDependencies/ClusterCommand.h"

#include <Rcpp.h>
//[[Rcpp::export]]
std::string MatrixToOpiMatrixCluster(const std::vector<int> &xPosition,
                        const std::vector<int> &yPosition, const std::vector<double> &data, const double cutoff,
                        const int iterations = 2)
{
    OptimatrixAdapter adapter(cutoff);
    const auto optiMatrix = adapter.ConvertToOptimatrix(xPosition,yPosition,data);
    auto* command = new ClusterCommand();
    command->SetMaxIterations(iterations);
    return command->runOptiCluster(optiMatrix);
}
#endif
int main() {
    // OptimatrixAdapterTestFixture fixture;
    // const auto xVals = std::vector<int>{0,0,0,1,1,2,3};
    // const auto yVals = std::vector<int>{1,2,4,2,4,4,4};
    // const auto data = std::vector<double>{0.02,0.04,0.025,0.01,0.028,0.045,0.05};
    // bool result = fixture.TestOptimatrixClosenessReturnsCorrectValue(xVals, yVals, data, 4);
    return 0;
}