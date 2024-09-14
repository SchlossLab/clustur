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
#include "MothurDependencies/ColumnDistanceMatrixReader.h"
#include "MothurDependencies/SharedFileBuilder.h"

#if DEBUG_RCPP
#include <Rcpp.h>

//[[Rcpp::export]]
void WritePhylipFile(const std::vector<int> &xPosition,
                      const std::vector<int> &yPosition, const std::vector<double> &data,
                      const double cutoff, const Rcpp::DataFrame& countTable, const std::string& saveLocation) {

    CountTableAdapter countTableAdapter;
    countTableAdapter.CreateDataFrameMap(countTable);
    MatrixAdapter adapter(xPosition, yPosition, data, cutoff, false, countTableAdapter);
    adapter.CreatePhylipFile(saveLocation);
}

//[[Rcpp::export]]
void WriteColumnFile(const std::vector<int> &xPosition,
                      const std::vector<int> &yPosition, const std::vector<double> &data,
                      const double cutoff, const Rcpp::DataFrame& countTable, const std::string& saveLocation) {

    CountTableAdapter countTableAdapter;
    countTableAdapter.CreateDataFrameMap(countTable);
    MatrixAdapter adapter(xPosition, yPosition, data, cutoff, false, countTableAdapter);
    adapter.CreateColumnDataFile(saveLocation, cutoff);
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
                                                  const int maxIterations = 100, const bool shuffle = true,
                                                  const bool isSim = false) {
    CountTableAdapter countTableAdapter;
    countTableAdapter.CreateDataFrameMap(countTable);
    MatrixAdapter adapter(xPosition, yPosition, data, cutoff, isSim, countTableAdapter);
    const auto sparseMatix = adapter.DistanceMatrixToSquareMatrix();
    OptimatrixAdapter optiAdapter(cutoff);
    const auto optiMatrix = optiAdapter.ConvertToOptimatrix(sparseMatix, isSim);
    ClusterCommand command;

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
                           const Rcpp::DataFrame& countTable,
                           const bool isSimularity) {
    CountTableAdapter countTableAdapter;
    countTableAdapter.CreateDataFrameMap(countTable);
    MatrixAdapter adapter(xPosition, yPosition, data, cutoff, isSimularity, countTableAdapter);
    ClusterCommand command;
    const auto sparseMatix = adapter.CreateSparseMatrix();
    const auto listVector = adapter.GetListVector();

    const auto result = command.runMothurCluster(method, sparseMatix, cutoff, listVector);
    Rcpp::DataFrame clusterDataFrame = result->GetListVector().listVector->CreateDataFrameFromList(result->GetListVector().label);
    Rcpp::DataFrame tidySharedDataFrame = CreateSharedDataFrame(countTableAdapter, result);
    delete(result);
    delete(listVector);
    return {tidySharedDataFrame, clusterDataFrame};

}

//[[Rcpp::export]]
std::vector<Rcpp::DataFrame> ClusterWithPhylip(const std::string& phylipFilePath,
                           const double cutoff,
                           const std::string& method,
                           const Rcpp::DataFrame& countTable,
                           const bool isSimularity) {
    CountTableAdapter countTableAdapter;
    countTableAdapter.CreateDataFrameMap(countTable);
    ReadPhylipMatrix reader(cutoff, isSimularity);
    reader.read(phylipFilePath);
    ClusterCommand command;
    const auto sparseMatix = reader.getDMatrix();
    const auto listVector = reader.getListVector();
    const auto result = command.runMothurCluster(method, sparseMatix, cutoff, listVector);
    Rcpp::DataFrame clusterDataFrame = result->GetListVector().listVector->CreateDataFrameFromList(result->GetListVector().label);
    Rcpp::DataFrame tidySharedDataFrame = CreateSharedDataFrame(countTableAdapter, result);
    delete(result);
    delete(listVector);
    return {tidySharedDataFrame, clusterDataFrame};
}

//[[Rcpp::export]]
std::vector<Rcpp::DataFrame> ClusterWithColumn(const std::string& columnFilePath,
                           const double cutoff,
                           const std::string& method,
                           const Rcpp::DataFrame& countTable,
                           const bool isSimularity) {
    CountTableAdapter countTableAdapter;
    countTableAdapter.CreateDataFrameMap(countTable);
    ColumnDistanceMatrixReader colReader(cutoff, isSimularity);
    colReader.Read(countTableAdapter, columnFilePath);
    ClusterCommand command;
    const auto sparseMatix = colReader.GetSparseDataMatrix();
    const auto listVector = colReader.GetListVector();
    const auto result = command.runMothurCluster(method, sparseMatix, cutoff, listVector);
    Rcpp::DataFrame clusterDataFrame = result->GetListVector().listVector->CreateDataFrameFromList(result->GetListVector().label);
    Rcpp::DataFrame tidySharedDataFrame = CreateSharedDataFrame(countTableAdapter, result);
    delete(result);
    delete(listVector);
    return {tidySharedDataFrame, clusterDataFrame};
}


//[[Rcpp::export]]
std::vector<Rcpp::DataFrame> OptiClusterPhylip(const std::string& filePath,
                                                const double cutoff, const Rcpp::DataFrame& countTable,
                                                const int maxIterations = 100, const bool shuffle = true,
                                                const bool isSim = false) {
    CountTableAdapter countTableAdapter;
    countTableAdapter.CreateDataFrameMap(countTable);
    ReadPhylipMatrix reader;
    const auto sparseMatix = reader.readToRowData(filePath);
    OptimatrixAdapter optiAdapter(cutoff);
    const auto optiMatrix = optiAdapter.ConvertToOptimatrix(sparseMatix, isSim);
    ClusterCommand command;

    command.SetOpticlusterRandomShuffle(shuffle);
    command.SetMaxIterations(maxIterations);
    const auto* result = command.runOptiCluster(optiMatrix, cutoff);
    Rcpp::DataFrame clusterDataFrame = result->GetListVector().listVector->CreateDataFrameFromList(result->GetListVector().label);
    Rcpp::DataFrame tidySharedDataFrame = CreateSharedDataFrame(countTableAdapter, result);
    delete(result);
    return {tidySharedDataFrame, command.GetSensitivityData(), command.GetClusterMetrics(), clusterDataFrame};
}

//[[Rcpp::export]]
std::vector<Rcpp::DataFrame> OptiClusterColumnDist(const std::string& filePath,
                                                const double cutoff, const Rcpp::DataFrame& countTable,
                                                const int maxIterations = 100, const bool shuffle = true,
                                                const bool isSim = false) {
    CountTableAdapter countTableAdapter;
    countTableAdapter.CreateDataFrameMap(countTable);
    ColumnDistanceMatrixReader reader(cutoff, isSim);
    const auto sparseMatix = reader.readToRowData(countTableAdapter, filePath);
    OptimatrixAdapter optiAdapter(cutoff);
    const auto optiMatrix = optiAdapter.ConvertToOptimatrix(sparseMatix, isSim);
    ClusterCommand command;

    command.SetOpticlusterRandomShuffle(shuffle);
    command.SetMaxIterations(maxIterations);
    const auto* result = command.runOptiCluster(optiMatrix, cutoff);
    Rcpp::DataFrame clusterDataFrame = result->GetListVector().listVector->CreateDataFrameFromList(result->GetListVector().label);
    Rcpp::DataFrame tidySharedDataFrame = CreateSharedDataFrame(countTableAdapter, result);
    delete(result);
    return {tidySharedDataFrame, command.GetSensitivityData(), command.GetClusterMetrics(), clusterDataFrame};
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
