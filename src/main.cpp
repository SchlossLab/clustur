#include <iostream>
#include <fstream>
#include <vector>
//TODO Next week: We have to separate the code out into a package! We have all the source code, now need to consider how
//it will look in a package!
#include "Adapters/OptimatrixAdapter.h"
#include "TestHelpers/TestHelper.h"
#include "Adapters/MatrixAdapter.h"
#include "MothurDependencies/ClusterCommand.h"
#include "MothurDependencies/RSparseMatrix.h"
#include "MothurDependencies/OptiMatrix.h"

#include "Adapters/CountTableAdapter.h"
#include "MothurDependencies/SharedFileBuilder.h"
#include "Tests/DataFrameAdapterTestFixture.h"

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
    SharedFileBuilder builder;
    std::unordered_map<std::string, RAbundVector> map;
    std::unordered_map<std::string, ListVector> listMap;
    const auto listVectors = result->GetListVector();
    listMap[listVectors.label] = *listVectors.listVector;
    const SharedFile* sharedFile = builder.BuildSharedFile(listMap, countTable);
    Rcpp::DataFrame tidySharedDataFrame = sharedFile->PrintData();
    delete(sharedFile);
    return tidySharedDataFrame;
}
//[[Rcpp::export]]
std::vector<Rcpp::DataFrame> MatrixToOpiMatrixCluster(const std::vector<int> &xPosition,
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
    const auto* result = command.runOptiCluster(optiMatrix, cutoff);
    Rcpp::DataFrame clusterDataFrame = result->GetListVector().listVector->CreateDataFrameFromList(result->GetListVector().label);
    Rcpp::DataFrame tidySharedDataFrame = CreateSharedDataFrame(countTableAdapter, result);
    delete(result);
    return {tidySharedDataFrame, command.GetSensitivityData(), command.GetClusterMetrics(), clusterDataFrame};
}


//[[Rcpp::export]]
std::vector<Rcpp::DataFrame> ClassicCluster(const std::vector<int> &xPosition,
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
    Rcpp::DataFrame clusterDataFrame = result->GetListVector().listVector->CreateDataFrameFromList(result->GetListVector().label);
    Rcpp::DataFrame tidySharedDataFrame = CreateSharedDataFrame(countTableAdapter, result);
    delete(result);
    delete(listVector);
    return {tidySharedDataFrame, clusterDataFrame};

}

#endif
// int main() {
//     CountTableAdapterTestFixture testFixture;
//     bool res = testFixture.TestGetGroups({"no-group"});
//     res = testFixture.TestGetGroups({});
//     res = testFixture.TestGetGroups({"1"});
//
//
//
//     return 0;
// }
