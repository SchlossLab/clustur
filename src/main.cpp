#include <iostream>
#include <fstream>
//TODO Next week: We have to separate the code out into a package! We have all the source code, now need to consider how
//it will look in a package!
#include "Adapters/OptimatrixAdapter.h"
#include "TestHelpers/TestHelper.h"
#include "Adapters/MatrixAdapter.h"
#include "MothurDependencies/ClusterCommand.h"
#include "MothurDependencies/CompleteLinkage.h"
#include "Tests/PhylipReaderTestFixture.h"

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
    PhylipReaderTestFixture fixture;
    MatrixAdapter adapter({1,2,3,4,5}, {2,3,4,5,6}, {.1,.11,.12,.15,.25}, 0.2);
    fixture.TestGetListVector(adapter.DistanceMatrixToSquareMatrix(), true);
}
