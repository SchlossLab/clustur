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

#include "Adapters/CountTableAdapter.h"
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

Rcpp::DataFrame CreateSharedDataFrame(const CountTableAdapter& countTable, const ClusterExport* result) {
    std::string clusterResults = result->Print();
    SharedFileBuilder builder;
    std::unordered_map<std::string, RAbundVector> map;
    std::unordered_map<std::string, ListVector> listMap;
    auto listVectors = result->GetListVectors();
    for(int i = 0; i < listVectors.size(); i++) {
        ListVector vector = *listVectors[i];
        listMap[result->GetLabel(i)] = vector;
    }
    SharedFile* sharedFile = builder.BuildSharedFile(listMap, countTable);
    Rcpp::DataFrame tidySharedDataFrame = sharedFile->PrintData();

    delete(sharedFile);
    return tidySharedDataFrame;
}
//[[Rcpp::export]]
Rcpp::DataFrame MatrixToOpiMatrixCluster(const std::vector<int> &xPosition,
                                                  const std::vector<int> &yPosition, const std::vector<double> &data,
                                                  const double cutoff, const Rcpp::DataFrame& countTable,
                                                  const int maxIterations = 100, const bool shuffle = true) {
    OptimatrixAdapter adapter(cutoff);
    const auto optiMatrix = adapter.ConvertToOptimatrix(xPosition, yPosition, data);
    ClusterCommand command;
    CountTableAdapter countTableAdapter;
    countTableAdapter.CreateDataFrameMap(countTable);
    command.SetOpticlusterRandomShuffle(shuffle);
    command.SetMaxIterations(maxIterations);
    const auto* result = command.runOptiCluster(optiMatrix);
    Rcpp::DataFrame tidySharedDataFrame = CreateSharedDataFrame(countTableAdapter, result);
    delete(result);
    return tidySharedDataFrame;
    //return {clusterResults, command.GetSensitivityData(), command.GetClusterMetrics()} ;
}


//[[Rcpp::export]]
Rcpp::DataFrame ClassicCluster(const std::vector<int> &xPosition,
                           const std::vector<int> &yPosition,
                           const std::vector<double> &data,
                           const double cutoff,
                           const std::string& method,
                           const Rcpp::DataFrame& countTable) {
    MatrixAdapter adapter(xPosition, yPosition, data, cutoff);
    ClusterCommand command;
    const auto sparseMatix = adapter.CreateSparseMatrix();
    const auto listVector = adapter.GetListVector();
    CountTableAdapter countTableAdapter;
    countTableAdapter.CreateDataFrameMap(countTable);
    const auto result = command.runMothurCluster(method, sparseMatix, cutoff, listVector);
    Rcpp::DataFrame tidySharedDataFrame = CreateSharedDataFrame(countTableAdapter, result);

    delete(result);
    delete(listVector);
    return tidySharedDataFrame;

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
