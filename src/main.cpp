#include <iostream>
#include <fstream>
#include <vector>
//TODO Next week: We have to separate the code out into a package! We have all the source code, now need to consider how
//it will look in a package!
#include "Adapters/OptimatrixAdapter.h"
#include "TestHelpers/TestHelper.h"
#include "Adapters/MatrixAdapter.h"
#include "MothurDependencies/ClusterCommand.h"
#include "MothurDependencies/OptiMatrix.h"
#include "Adapters/CountTableAdapter.h"
#include "Adapters/DistanceFilesToSparse.h"
#include "MothurDependencies/ColumnDistanceMatrixReader.h"
#include "MothurDependencies/SharedFileBuilder.h"
#include "Adapters/DistanceFileReader.h"
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
                                                  const bool isSim = false) {
    CountTableAdapter countTableAdapter;
    countTableAdapter.CreateDataFrameMap(countTable);
    MatrixAdapter adapter(xPosition, yPosition, data, cutoff, isSim, countTableAdapter);
    const auto sparseMatix = adapter.DistanceMatrixToSquareMatrix();
    OptimatrixAdapter optiAdapter(cutoff);
    const auto optiMatrix = optiAdapter.ConvertToOptimatrix(sparseMatix, isSim);
    ClusterCommand command;
    const auto* result = command.runOptiCluster(optiMatrix, cutoff);
    Rcpp::DataFrame clusterDataFrame = result->GetListVector().listVector->CreateDataFrameFromList(result->GetListVector().label);
    Rcpp::DataFrame tidySharedDataFrame = CreateSharedDataFrame(countTableAdapter, result);
    delete(result);
    return {tidySharedDataFrame, clusterDataFrame,command.GetSensitivityData(), command.GetClusterMetrics()};
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
    reader.Read(phylipFilePath);
    ClusterCommand command;
    const auto sparseMatix = reader.GetSparseMatrix();
    const auto listVector = reader.GetListVector();
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
    const auto sparseMatix = colReader.GetSparseMatrix();
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
                                                const bool isSim = false) {
    CountTableAdapter countTableAdapter;
    countTableAdapter.CreateDataFrameMap(countTable);
    ReadPhylipMatrix reader;
    const auto sparseMatix = reader.ReadToRowData(filePath);
    OptimatrixAdapter optiAdapter(cutoff);
    const auto optiMatrix = optiAdapter.ConvertToOptimatrix(sparseMatix, isSim);
    ClusterCommand command;
    const auto* result = command.runOptiCluster(optiMatrix, cutoff);
    Rcpp::DataFrame clusterDataFrame = result->GetListVector().listVector->CreateDataFrameFromList(result->GetListVector().label);
    Rcpp::DataFrame tidySharedDataFrame = CreateSharedDataFrame(countTableAdapter, result);
    delete(result);
    return {tidySharedDataFrame, clusterDataFrame,command.GetSensitivityData(), command.GetClusterMetrics()};
}

//[[Rcpp::export]]
std::vector<Rcpp::DataFrame> OptiClusterColumnDist(const std::string& filePath,
                                                const double cutoff, const Rcpp::DataFrame& countTable,
                                                const bool isSim = false) {
    CountTableAdapter countTableAdapter;
    countTableAdapter.CreateDataFrameMap(countTable);
    ColumnDistanceMatrixReader reader(cutoff, isSim);
    const auto sparseMatix = reader.ReadToRowData(countTableAdapter, filePath);
    OptimatrixAdapter optiAdapter(cutoff);
    const auto optiMatrix = optiAdapter.ConvertToOptimatrix(sparseMatix, isSim);
    ClusterCommand command;
    const auto* result = command.runOptiCluster(optiMatrix, cutoff);
    Rcpp::DataFrame clusterDataFrame = result->GetListVector().listVector->CreateDataFrameFromList(result->GetListVector().label);
    Rcpp::DataFrame tidySharedDataFrame = CreateSharedDataFrame(countTableAdapter, result);
    delete(result);
    return {tidySharedDataFrame, clusterDataFrame,command.GetSensitivityData(), command.GetClusterMetrics()};
}

//[[Rcpp::export]]
std::vector<std::vector<double>> DistanceFileToSparseMatrix(const Rcpp::DataFrame& countTable, const std::string& filePath,
    const std::string& method) {
    CountTableAdapter countTableAdapter;
    countTableAdapter.CreateDataFrameMap(countTable);
    const DistanceFilesToSparse sparseConverter(countTableAdapter, method, filePath);
    return sparseConverter.CreateSparseList();
}

//[[Rcpp::export]]
SEXP ProcessDistanceFiles(const std::string& filePath, const Rcpp::DataFrame& countTable, double cutoff, bool isSim) {
    std::fstream data(filePath);
    std::unordered_map<bool, std::string> map;
    map[true] = "This is a phylip file. Processing now...";
    map[false] = "This is a column file. Processing now...";

    if(!data.is_open()) {
        Rcpp::Rcout << "Please enter a valid file path\n";
        Rcpp::stop("Invalid file path.");
    }
    std::string line;
    std::getline(data, line);
    bool isPhylip = true;
    int count = 0;
    for(const auto& character: line) {
        if(character == '\t')
            count++;
    }
    if(count > 1)
        isPhylip = false;
    Rcpp::Rcout << map[isPhylip] << "\n";
    data.close();

    CountTableAdapter adapter;
    adapter.CreateDataFrameMap(countTable);
    if(isPhylip) {
        DistanceFileReader* read = new ReadPhylipMatrix(cutoff, isSim);
        read->SetCountTable(adapter);
        read->SetSparseMatrix(read->ReadToRowData(filePath));
        return Rcpp::XPtr<DistanceFileReader>(read);
    }
    DistanceFileReader* read = new ColumnDistanceMatrixReader(cutoff, isSim);
    read->SetCountTable(adapter);
    read->SetSparseMatrix(read->ReadToRowData(adapter, filePath));
    return Rcpp::XPtr<DistanceFileReader>(read);
}

SEXP ReadSparseMatrix(const CountTableAdapter& countTableAdapter, const std::vector<int> &xPosition,
    const std::vector<int> &yPosition, const std::vector<double> &data,
    const double cutoff, const bool isSim) {
    MatrixAdapter adapter(xPosition, yPosition, data, cutoff, isSim, countTableAdapter);
    const auto sparseMatix = adapter.CreateSparseMatrix();
    const auto listVector = adapter.GetListVector();
}

//[[Rcpp::export]]
Rcpp::DataFrame GetDistanceDataFrame(const SEXP& fileReader) {
    const Rcpp::XPtr<DistanceFileReader> ptr(fileReader);
    return ptr.get()->SparseMatrixToDataFrame();
}

//[[Rcpp::export]]
Rcpp::DataFrame GetCountTable(const SEXP& fileReader) {
    const Rcpp::XPtr<DistanceFileReader> ptr(fileReader);
    return ptr.get()->GetCountTable();
}
#endif
