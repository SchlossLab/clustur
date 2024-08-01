#include <iostream>
#include <fstream>
//TODO Next week: We have to separate the code out into a package! We have all the source code, now need to consider how
//it will look in a package!
#include "Adapters/OptimatrixAdapter.h"
#include "TestHelpers/TestHelper.h"
#include "Adapters/MatrixAdapter.h"
#include "MothurDependencies/ClusterCommand.h"
#include "MothurDependencies/RSparseMatrix.h"
#include "MothurDependencies/OptiMatrix.h"
#include <vector>

#include "MothurDependencies/SharedFileBuilder.h"

#if DEBUG_RCPP
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
    SharedFileBuilder builder(RSpraseMatrix(xPosition, yPosition, data));
    builder.CreateNameMap();
    builder.CreateCountTable();
    builder.OutputFiles();
    ClusterCommand command;
    //Race Condition, going to have to look for a fix in the future, but Create Sparse Matrix has to go first
    const auto sparseMatix = adapter.CreateSparseMatrix();
    const auto listVector = adapter.GetListVector();
    const auto result = command.runMothurCluster(method, sparseMatix, cutoff, listVector);
    const auto shared = builder.BuildSharedFile(listVector);
    shared->PrintData();
    delete(shared);
    return result;

}
#endif
// int main() {
//     const auto xVals = std::vector<int>{1,1,1,2,2,3,5};
//     const auto yVals = std::vector<int>{2,3,5,3,5,5,5};
//     const auto data = std::vector<double>{0.02,0.04,0.25,0.1,0.28,0.05,0.15};
//     ClassicCluster(xVals,yVals,data,0.2, "furthest");
//     // SharedFileBuilder builder(RSpraseMatrix(xVals, yVals, data));
//     // builder.CreateNameMap();
//     // builder.CreateCountTable();
//     // builder.OutputFiles();
//     return 0;
// }