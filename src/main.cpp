#include <iostream>
#include <fstream>
//TODO Next week: We have to separate the code out into a package! We have all the source code, now need to consider how
//it will look in a package!
#include "MothurDependencies/Cluster.h"
#include "MothurDependencies/SingleLinkage.h"
#include "TestHelpers/TestHelper.h"
#include "Tests/UtilsTestFixture.h"
#if DEBUG_RCPP
#include <vector>
#include "Adapters/OptimatrixAdapter.h"
#include "MothurDependencies/OptiMatrix.h"
#include "MothurDependencies/ClusterCommand.h"

#include <Rcpp.h>
//[[Rcpp::export]]
std::vector<std::string> MatrixToOpiMatrixCluster(const std::vector<int> &xPosition,
                        const std::vector<int> &yPosition, const std::vector<double> &data, const double cutoff,
                        const int iterations = 2, const bool shuffle = true)
{
    OptimatrixAdapter adapter(cutoff);
    const auto optiMatrix = adapter.ConvertToOptimatrix(xPosition,yPosition,data);
    auto* command = new ClusterCommand();
    command ->SetOpticlusterRandomShuffle(shuffle);
    command->SetMaxIterations(iterations);
    return command->runOptiCluster(optiMatrix);
}
#endif
SparseDistanceMatrix* MatrixToCluster(const std::vector<int> &xPosition,
                        const std::vector<int> &yPosition, const std::vector<double> &data, const double cutoff) {
    const size_t count = data.size();
    std::vector<PDistCell> cells(count);
    auto* distanceMatrix = new SparseDistanceMatrix();
    distanceMatrix->resize(count);
    for(size_t i = 0; i < data.size(); i++) {

        cells[i] = PDistCell(xPosition[i] , data[i]);
        distanceMatrix->addCell(i, PDistCell(xPosition[i] - 1, data[i]));
    }
    return distanceMatrix;
}

int main() {
    double cutoff = 0.2;
    auto distanceMatrix = MatrixToCluster({1,2,3,4}, {2,3,4,5}, {5,10,15,20}, cutoff);
    SingleLinkage* cluster = new SingleLinkage(new RAbundVector(), new ListVector(),
        distanceMatrix, 0.2, "nearest", cutoff);
    cluster->update(cutoff);
    return 0;
}