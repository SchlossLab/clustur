#include <iostream>
#include <fstream>
#include <vector>
#include "Adapters/OptimatrixAdapter.h"
#include "Adapters/MatrixAdapter.h"
#include "MothurDependencies/ClusterCommand.h"
#include "MothurDependencies/ListVector.h"
#include "MothurDependencies/OptiMatrix.h"
#include "Adapters/CountTableAdapter.h"
#include "MothurDependencies/ColumnDistanceMatrixReader.h"
#include "MothurDependencies/SharedFileBuilder.h"
#include "Adapters/DistanceFileReader.h"
#include <Rcpp.h>
#include <cctype>


Rcpp::DataFrame CreateSharedDataFrame(const CountTableAdapter& countTable, const ClusterExport* result,
    const std::string& binName) {
    SharedFileBuilder builder;
    std::unordered_map<std::string, RAbundVector> map;
    std::unordered_map<std::string, ListVector> listMap;
    const ListVectorPair listVectors = result->GetListVector();
    const SharedFile* sharedFile = builder.BuildSharedFile(*listVectors.listVector, countTable, binName);
    Rcpp::DataFrame tidySharedDataFrame = sharedFile->PrintData(binName);
    delete(sharedFile);
    return tidySharedDataFrame;
}

//[[Rcpp::export]]
bool DetermineIfPhylipOrColumnFile(const std::string& filePath) {
    std::ifstream data(filePath);

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
    if(isPhylip) {
        DistanceFileReader* read = new ReadPhylipMatrix(cutoff, isSim);
        read->CreateCountTableAdapter(countTable);
        read->Read(filePath);
        return Rcpp::XPtr<DistanceFileReader>(read);
    }
    DistanceFileReader* read = new ColumnDistanceMatrixReader(cutoff, isSim);
    read->CreateCountTableAdapter(countTable);
    read->Read(filePath);
    return Rcpp::XPtr<DistanceFileReader>(read);
}

//[[Rcpp::export]]
SEXP ProcessSparseMatrix(const std::vector<int> &xPosition,
    const std::vector<int> &yPosition, const std::vector<double> &data, const Rcpp::DataFrame& countTable,
    const double cutoff, const bool isSim) {
    CountTableAdapter countTableAdapter;
    countTableAdapter.CreateDataFrameMap(countTable);
    MatrixAdapter adapter(xPosition, yPosition, data, cutoff, isSim, countTableAdapter);
    auto* sparseDistanceMatrix = new SparseDistanceMatrix(adapter.CreateSparseMatrix());
    auto* listVec =  new ListVector(adapter.CreateListVector());
    auto* read = new DistanceFileReader(sparseDistanceMatrix,listVec,cutoff, isSim);
    read->CreateCountTableAdapter(countTable);
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
Rcpp::List Cluster(const SEXP& DistanceData,const std::string& method, const std::string& featureColumnName,
    const std::string& binColumnName, const double cutoff) {
    const Rcpp::XPtr<DistanceFileReader> distanceData(DistanceData);
    const CountTableAdapter countTableAdapter = distanceData.get()->GetCountTableAdapter();
    ClusterCommand command;
    const auto lastCutoff = distanceData.get()->GetCutoff();
    const auto listVector = distanceData.get()->GetListVector(); // Going to have to make a copy of list vector, this two values are definitely being changed
    auto sparseMatrix = distanceData.get()->GetSparseMatrix(); // Going to have to make a copy of sparse matrix
    if(cutoff < lastCutoff)
        sparseMatrix->FilterSparseMatrix(cutoff);
    const auto result = command.runMothurCluster(method, sparseMatrix, cutoff, listVector);
    const auto label = result->GetListVector().label;
    const Rcpp::DataFrame clusterDataFrame = result->GetListVector().listVector->CreateDataFrameFromList(
        featureColumnName, binColumnName);
    const Rcpp::DataFrame tidySharedDataFrame = CreateSharedDataFrame(countTableAdapter, result, binColumnName);
    delete(result);
    delete(listVector);
    delete(sparseMatrix);
    return Rcpp::List::create(Rcpp::Named("label") = std::stod(label),
    Rcpp::Named("abundance") = tidySharedDataFrame,
    Rcpp::Named("cluster") = clusterDataFrame);
}

//[[Rcpp::export]]
Rcpp::List OptiCluster(const SEXP& DistanceData, const std::string& featureColumnName, const std::string& binColumnName,
    const double cutoff) {
    const Rcpp::XPtr<DistanceFileReader> distanceData(DistanceData);
    const CountTableAdapter countTableAdapter = distanceData.get()->GetCountTableAdapter();
    const auto sparseMatix =  distanceData.get()->GetSparseMatrix();
    const auto listVector = distanceData.get()->GetListVector();
    const bool isSim = distanceData.get()->GetIsSimularity();
    OptimatrixAdapter optiAdapter(cutoff);
    const auto optiMatrix = optiAdapter.ConvertToOptimatrix(sparseMatix, listVector, isSim);
    delete(sparseMatix);
    delete(listVector);
    ClusterCommand command;
    const auto* result = command.runOptiCluster(optiMatrix, cutoff);
    const auto label = result->GetListVector().label;
    const Rcpp::DataFrame clusterDataFrame = result->GetListVector().listVector->CreateDataFrameFromList(
        featureColumnName, binColumnName);
    const Rcpp::DataFrame tidySharedDataFrame = CreateSharedDataFrame(countTableAdapter, result, binColumnName);
    delete(result);
    return Rcpp::List::create(Rcpp::Named("label") = std::stod(label),
      Rcpp::Named("abundance") = tidySharedDataFrame,
      Rcpp::Named("cluster") = clusterDataFrame,
      Rcpp::Named("cluster_metrics") = command.GetSensitivityData(),
      Rcpp::Named("iteration_metrics") = command.GetClusterMetrics());
}

//[[Rcpp::export]]
Rcpp::DataFrame CreateDataFrameFromSparseCountTable(const Rcpp::DataFrame& countTable) {
    CountTableAdapter adapter;
    adapter.CreateDataFrameMapFromSparseCountTable(countTable);
    return adapter.ReCreateDataFrame();
}
