#include <iostream>
#include <fstream>
//TODO Next week: We have to separate the code out into a package! We have all the source code, now need to consider how
//it will look in a package!
#include "Adapters/OptimatrixAdapter.h"
#include "TestHelpers/TestHelper.h"
#include "Adapters/MatrixAdapter.h"
#include "MothurDependencies/ClusterCommand.h"
#include "MothurDependencies/CompleteLinkage.h"
#include "Tests/MatrixAdapterTestFixture.h"
#include "Tests/PhylipReaderTestFixture.h"
#include "Tests/RAbundVectorTestFixture.h"

#if DEBUG_RCPP
#include <vector>
#include "Adapters/OptimatrixAdapter.h"
#include "MothurDependencies/OptiMatrix.h"
#include "MothurDependencies/ClusterCommand.h"
#include <Rcpp.h>

//[[Rcpp::export]]
void WritePhylipFile(const std::vector<int> &xPosition,
                      const std::vector<int> &yPosition, const std::vector<double> &data,
                      const double cutoff, const std::string& saveLocation) {
    MatrixAdapter adapter(xPosition, yPosition, data, cutoff);
    adapter.CreatePhylipFile(saveLocation);
}

//[[Rcpp::export]]
std::vector<std::string> MatrixToOpiMatrixCluster(const std::vector<int> &xPosition,
                                                  const std::vector<int> &yPosition, const std::vector<double> &data,
                                                  const double cutoff,
                                                  const int maxIterations = 100, const bool shuffle = true) {
    OptimatrixAdapter adapter(cutoff);
    const auto optiMatrix = adapter.ConvertToOptimatrix(xPosition, yPosition, data);
    ClusterCommand command;
    command.SetOpticlusterRandomShuffle(shuffle);
    command.SetMaxIterations(maxIterations);
    return command.runOptiCluster(optiMatrix);
}


//[[Rcpp::export]]
std::string ClassicCluster(const std::vector<int> &xPosition,
                           const std::vector<int> &yPosition,
                           const std::vector<double> &data,
                           const double cutoff,
                           const std::string& method) {
    MatrixAdapter adapter(xPosition, yPosition, data, cutoff);
    ClusterCommand command;
    return command.runMothurCluster(method, adapter.CreateSparseMatrix(), cutoff, adapter.GetListVector());}
#endif

int main() {
    RAbundVectorTestFixture fixture;
    bool result = fixture.TestClear(true);
    result = fixture.TestGet(1, 1);
    result = fixture.TestQuicksort(true);
    result = fixture.TestRemove(1,true);
    result = fixture.TestResize(3,true);
    result = fixture.TestSet(1,2, true);
    result = fixture.TestSize(6);
    result = fixture.TestGetOverload({1,1,1,1,1,1});
    result = fixture.TestPushBack(2, true);
    result = fixture.TestGetMaxRank(1);
    result = fixture.TestGetNumBins(6);
    result = fixture.TestGetNumSeqs(6);
}
