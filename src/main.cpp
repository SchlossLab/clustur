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
#include "MothurDependencies/ColumnDistanceMatrixReader.h"
#include "MothurDependencies/SharedFileBuilder.h"
#include "Adapters/DistanceFileReader.h"
#if DEBUG_RCPP
#include <Rcpp.h>
#include <cctype>

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
    const ListVectorPair listVectors = result->GetListVector();
    const SharedFile* sharedFile = builder.BuildSharedFile(*listVectors.listVector, countTable);
    Rcpp::DataFrame tidySharedDataFrame = sharedFile->PrintData();
    delete(sharedFile);
    return tidySharedDataFrame;
}


//[[Rcpp::export]]
bool DetermineIfPhylipOrColumnFile(const std::string& filePath) {
    std::fstream data(filePath);

    if(!data.is_open()) {
        Rcpp::Rcout << "Please enter a valid file path\n";
        Rcpp::stop("Invalid file path.");
    }
    std::string line;
    std::getline(data, line);
    bool isPhylip = true;
    int count = 0;
    for(const auto& character: line) {
        if(character != '\n' && std::isspace(character))
            count++;
    }
    if(count > 1)
        isPhylip = false;
    data.close();
    return isPhylip;
}

//[[Rcpp::export]]
SEXP ProcessDistanceFiles(const std::string& filePath, const Rcpp::DataFrame& countTable, const double cutoff,
    const bool isSim) {
    const bool isPhylip = DetermineIfPhylipOrColumnFile(filePath);

    CountTableAdapter adapter;
    adapter.CreateDataFrameMap(countTable);
    if(isPhylip) {
        DistanceFileReader* read = new ReadPhylipMatrix(cutoff, isSim);
        const std::vector<RowData> rowDataMatrix = read->ReadToRowData(filePath);
        read->SetCountTable(adapter);
        read->SetRowDataMatrix(rowDataMatrix);
        read->ReadRowDataMatrix(rowDataMatrix);
        return Rcpp::XPtr<DistanceFileReader>(read);
    }
    DistanceFileReader* read = new ColumnDistanceMatrixReader(cutoff, isSim);
    const std::vector<RowData> rowDataMatrix = read->ReadToRowData(adapter, filePath);
    read->SetCountTable(adapter);
    read->SetRowDataMatrix(rowDataMatrix);
    read->ReadRowDataMatrix(rowDataMatrix);
    return Rcpp::XPtr<DistanceFileReader>(read);
}

//[[Rcpp::export]]
SEXP ProcessSparseMatrix(const std::vector<int> &xPosition,
    const std::vector<int> &yPosition, const std::vector<double> &data, const Rcpp::DataFrame& countTable,
    const double cutoff, const bool isSim) {
    CountTableAdapter countTableAdapter;
    countTableAdapter.CreateDataFrameMap(countTable);
    DistanceFileReader* read = new ReadPhylipMatrix(cutoff, isSim);
    MatrixAdapter adapter(xPosition, yPosition, data, cutoff, isSim, countTableAdapter);
    const std::vector<RowData> rowDataMatrix = adapter.DistanceMatrixToSquareMatrix();
    read->ReadRowDataMatrix(rowDataMatrix);
    read->SetRowDataMatrix(rowDataMatrix);
    read->SetCountTable(countTableAdapter);
    return Rcpp::XPtr<DistanceFileReader>(read);
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

//[[Rcpp::export]]
Rcpp::List Cluster(const SEXP& DistanceData, const std::string& method) {
    const Rcpp::XPtr<DistanceFileReader> distanceData(DistanceData);
    const CountTableAdapter countTableAdapter = distanceData.get()->GetCountTableAdapter();
    ClusterCommand command;
    const double cutoff = distanceData.get()->GetCutoff();
    const auto sparseMatix = distanceData.get()->GetSparseMatrix(); // Going to have to make a copy of sparse matrix
    const auto listVector = distanceData.get()->GetListVector(); // Going to have to make a copy of list vector, this two values are definitely being changed
    const auto result = command.runMothurCluster(method, sparseMatix, cutoff, listVector);
    const auto label = result->GetListVector().label;
    const Rcpp::DataFrame clusterDataFrame = result->GetListVector().listVector->CreateDataFrameFromList();
    const Rcpp::DataFrame tidySharedDataFrame = CreateSharedDataFrame(countTableAdapter, result);
    delete(result);
    delete(listVector);
    return Rcpp::List::create(Rcpp::Named("label") = std::stod(label),
    Rcpp::Named("abundance") = tidySharedDataFrame,
    Rcpp::Named("cluster") = clusterDataFrame);
}

//[[Rcpp::export]]
Rcpp::List OptiCluster(const SEXP& DistanceData) {
    const Rcpp::XPtr<DistanceFileReader> distanceData(DistanceData);
    const CountTableAdapter countTableAdapter = distanceData.get()->GetCountTableAdapter();
    const std::vector<RowData> sparseMatix =  distanceData.get()->GetRowDataMatrix();
    const double cutoff = distanceData.get()->GetCutoff();
    const bool isSim = distanceData.get()->GetIsSimularity();
    OptimatrixAdapter optiAdapter(cutoff);
    const auto optiMatrix = optiAdapter.ConvertToOptimatrix(sparseMatix, isSim);
    ClusterCommand command;
    const auto* result = command.runOptiCluster(optiMatrix, cutoff);
    const auto label = result->GetListVector().label;
    const Rcpp::DataFrame clusterDataFrame = result->GetListVector().listVector->CreateDataFrameFromList();
    const Rcpp::DataFrame tidySharedDataFrame = CreateSharedDataFrame(countTableAdapter, result);
    delete(result);
    return Rcpp::List::create(Rcpp::Named("label") = std::stod(label),
      Rcpp::Named("abundance") = tidySharedDataFrame,
      Rcpp::Named("cluster") = clusterDataFrame,
      Rcpp::Named("cluster_metrics") = command.GetSensitivityData(),
      Rcpp::Named("iteration_metrics") = command.GetClusterMetrics());
}

//[[Rcpp::export]]
Rcpp::DataFrame CreateDataFrameFromSparse(const Rcpp::DataFrame& countTable) {
    CountTableAdapter adapter;
    adapter.CreateDataFrameMapFromSparseCountTable(countTable);
    return adapter.ReCreateDataFrame();
}
#endif
