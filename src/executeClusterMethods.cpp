//
// Created by Gregory Johnson on 8/21/26.
//

#include <Rcpp.h>

#include "Adapters/OptimatrixAdapter.h"
#include "Clusters/AverageLinkage.h"
#include "Clusters/ClusterMethod.h"
#include "Clusters/CompleteLinkage.h"
#include "Clusters/OptiCluster.h"
#include "Clusters/Optifitcluster.h"
#include "Clusters/SingleLinkage.h"
#include "Clusters/WeightedLinkage.h"
#include "Clusters/Metrics/mcc.h"
#include "Clusters/Metrics/tptn.h"
#include "DataStructures/SplitMatrix.h"
#include "FileReaders/DistanceFileReader.h"
#include "SharedFileData/SharedFileBuilder.h"


Rcpp::DataFrame CreateSharedDataFrame(const CountTableAdapter& countTable, const ClusterExport* result,
                                      const std::string& binName) {
    std::unordered_map<std::string, RAbundVector> map;
    std::unordered_map<std::string, ListVector> listMap;
    const auto [listVector, label] = result->GetListVector();
    const SharedFile* sharedFile = SharedFileBuilder::BuildSharedFile(listVector, countTable, binName);
    Rcpp::DataFrame tidySharedDataFrame = sharedFile->PrintData(binName);
    delete(sharedFile);
    return tidySharedDataFrame;
}

//[[Rcpp::export]]
Rcpp::List Cluster(const SEXP& DistanceData,const std::string& method, const std::string& featureColumnName,
    const std::string& binColumnName, const double cutoff) {
    const Rcpp::XPtr<DistanceFileReader> distanceData(DistanceData);
    const CountTableAdapter countTableAdapter = distanceData.get()->GetCountTableAdapter();
    const auto lastCutoff = distanceData.get()->GetCutoff();
    ListVector* listVector = distanceData.get()->GetListVector(); // Going to have to make a copy of list vector, this two values are definitely being changed
    SparseDistanceMatrix* sparseMatrix = distanceData.get()->GetSparseMatrix(); // Going to have to make a copy of sparse matrix
    RAbundVector rAbund = listVector->getRAbundVector();
    if(cutoff < lastCutoff)
        sparseMatrix->FilterSparseMatrix(cutoff);
    ClusterMethod* clusterMethod = Utils::GetClusterMethod(method, listVector, sparseMatrix,
        rAbund, cutoff);
    const auto result = clusterMethod->Execute();
    const auto label = result->GetListVector().label;
    const Rcpp::DataFrame clusterDataFrame = result->GetListVector().listVector.CreateDataFrameFromList(
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
Rcpp::List OptiClust(const SEXP& DistanceData, const std::string& featureColumnName, const std::string& binColumnName,
    const double cutoff) {
    const Rcpp::XPtr<DistanceFileReader> distanceData(DistanceData);
    const CountTableAdapter countTableAdapter = distanceData.get()->GetCountTableAdapter();
    const auto sparseMatix =  distanceData.get()->GetSparseMatrix();
    const auto listVector = distanceData.get()->GetListVector();
    const bool isSim = distanceData.get()->GetIsSimularity();
    const OptimatrixAdapter optiAdapter(cutoff);
    auto* optiMatrix = optiAdapter.ConvertToOptimatrix(sparseMatix, listVector, isSim);
    delete(sparseMatix);
    delete(listVector);
    ClusterMetric* metric = nullptr;
    if (!optiMatrix->mccValidCalc()) {
        Rcpp::warning("[WARNING]: The mcc metric is not suitible for your data with a cutoff of " +
            std::to_string(cutoff) + " using tptn instead.");
        delete metric;
        metric = new TPTN();
    }
    else
        metric = new MCC();

    OptiCluster cluster(optiMatrix, metric, cutoff, 0);
    const auto* result = cluster.Execute();

    const Rcpp::DataFrame clusterMetricsDataFrame = cluster.GetSensitivityData();
    const Rcpp::DataFrame iterationsMetricsDataFrame = cluster.GetClusterMetrics();
    const auto label = result->GetListVector().label;
    const Rcpp::DataFrame clusterDataFrame = result->GetListVector().listVector.CreateDataFrameFromList(
        featureColumnName, binColumnName);
    const Rcpp::DataFrame tidySharedDataFrame = CreateSharedDataFrame(countTableAdapter, result, binColumnName);
    delete metric;
    delete result;
    delete optiMatrix;

    return Rcpp::List::create(Rcpp::Named("label") = std::stod(label),
      Rcpp::Named("abundance") = tidySharedDataFrame,
      Rcpp::Named("cluster") = clusterDataFrame,
      Rcpp::Named("cluster_metrics") = clusterMetricsDataFrame,
      Rcpp::Named("iteration_metrics") = iterationsMetricsDataFrame);
}


//[[Rcpp::export]]
Rcpp::List OptiFit(const SEXP& distData, const std::string& featureColumnName, const std::string& binColumnName,
    const double cutoff, const double fitPercent = 50, const bool isClosed = true, const bool printRef = true,
    const bool selfReference = false) {
    const Rcpp::XPtr<DistanceFileReader> distanceData(distData);
    const CountTableAdapter countTableAdapter = distanceData.get()->GetCountTableAdapter();
    const auto sparseMatix =  distanceData.get()->GetSparseMatrix();
    const auto listVector = distanceData.get()->GetListVector();
    const bool isSim = distanceData.get()->GetIsSimularity();
    const OptimatrixAdapter optiAdapter(cutoff);
    const auto* optiMatrix = optiAdapter.ConvertToOptimatrix(sparseMatix, listVector, isSim);
    auto* refMatrix = new OptiRefMatrix(optiMatrix, countTableAdapter, fitPercent, "");
    delete optiMatrix;
    delete(sparseMatix);
    delete(listVector);
    ClusterMetric* metric = new MCC();
    OptiFitCluster cluster(refMatrix, metric,"denovo", cutoff, 0, selfReference, printRef, isClosed);
    const auto* result = cluster.Execute();
    delete metric;
    delete refMatrix;
    const Rcpp::DataFrame clusterMetricsDataFrame = cluster.GetSensitivityData();
    const Rcpp::DataFrame iterationsMetricsDataFrame = cluster.GetClusterMetrics();
    const auto label = result->GetListVector().label;
    const Rcpp::DataFrame clusterDataFrame = result->GetListVector().listVector.CreateDataFrameFromList(
        featureColumnName, binColumnName);
    const Rcpp::DataFrame tidySharedDataFrame = CreateSharedDataFrame(countTableAdapter, result, binColumnName);
    delete(result);

    return Rcpp::List::create(Rcpp::Named("label") = std::stod(label),
      Rcpp::Named("abundance") = tidySharedDataFrame,
      Rcpp::Named("cluster") = clusterDataFrame,
      Rcpp::Named("cluster_metrics") = clusterMetricsDataFrame,
      Rcpp::Named("iteration_metrics") = iterationsMetricsDataFrame);
}

//[[Rcpp::export]]
Rcpp::List OptiFit2(const SEXP& distData, const std::string& featureColumnName, const std::string& binColumnName,
    const std::vector<std::string>& accnos,
    const double cutoff, const bool isClosed = true, const bool printRef = true, const bool selfReference = false) {
    const Rcpp::XPtr<DistanceFileReader> distanceData(distData);
    const CountTableAdapter countTableAdapter = distanceData.get()->GetCountTableAdapter();
    const auto sparseMatix =  distanceData.get()->GetSparseMatrix();
    const auto listVector = distanceData.get()->GetListVector();
    const bool isSim = distanceData.get()->GetIsSimularity();
    const OptimatrixAdapter optiAdapter(cutoff);
    const auto* optiMatrix = optiAdapter.ConvertToOptimatrix(sparseMatix, listVector, isSim);
    auto* refMatrix = new OptiRefMatrix(optiMatrix, countTableAdapter,
        {accnos.begin(), accnos.end()});
    delete optiMatrix;
    delete(sparseMatix);
    delete(listVector);
    ClusterMetric* metric = new MCC();
    OptiFitCluster cluster(refMatrix, metric,"userref", cutoff, 0, selfReference, printRef, isClosed);

    const auto* result = cluster.Execute();
    delete metric;
    delete refMatrix;
    const Rcpp::DataFrame clusterMetricsDataFrame = cluster.GetSensitivityData();
    const Rcpp::DataFrame iterationsMetricsDataFrame = cluster.GetClusterMetrics();
    const auto label = result->GetListVector().label;
    const Rcpp::DataFrame clusterDataFrame = result->GetListVector().listVector.CreateDataFrameFromList(
        featureColumnName, binColumnName);
    const Rcpp::DataFrame tidySharedDataFrame = CreateSharedDataFrame(countTableAdapter, result, binColumnName);
    delete(result);

    return Rcpp::List::create(Rcpp::Named("label") = std::stod(label),
      Rcpp::Named("abundance") = tidySharedDataFrame,
      Rcpp::Named("cluster") = clusterDataFrame,
      Rcpp::Named("cluster_metrics") = clusterMetricsDataFrame,
      Rcpp::Named("iteration_metrics") = iterationsMetricsDataFrame);
}

//[[Rcpp::export]]
Rcpp::List OptiFit3(const SEXP& combinedData, const Rcpp::DataFrame& refList, const std::vector<std::string>& accnos,
    const float fitPercent, const std::string& featureColumnName, const std::string& binColumnName,
    const double cutoff, const bool isClosed = true, const bool printRef = false, const bool selfReference = true) {
// fitPercent = fitPercent = ((count-refCount) / static_cast<float>(count));
    const ListVector refListOtuVector = Utils::CreateListVectorFromOtuList(refList["bin_name"],
        refList["sequence_name"]);
    const Rcpp::XPtr<DistanceFileReader> combinedDistanceData(combinedData);
    const CountTableAdapter combinedCountTableAdapter = combinedDistanceData.get()->GetCountTableAdapter();
    const SparseDistanceMatrix* combinedSparseMartix =  combinedDistanceData.get()->GetSparseMatrix();
    const ListVector* combinedListVector = combinedDistanceData.get()->GetListVector();
    const bool combinedIsSim = combinedDistanceData.get()->GetIsSimularity();
    const OptimatrixAdapter combinedOptiAdapter(cutoff);
    const auto* combinedOptiMatrix = combinedOptiAdapter.ConvertToOptimatrix(combinedSparseMartix, combinedListVector, combinedIsSim);
    delete combinedSparseMartix;
    delete combinedListVector;

    auto* refMatrix = new OptiRefMatrix(combinedOptiMatrix, combinedCountTableAdapter,
        {accnos.begin(), accnos.end()});
    ClusterMetric* metric = new MCC();
    OptiFitCluster cluster(refMatrix, metric, refListOtuVector, cutoff, 0, selfReference, printRef, isClosed);
    const auto* result = cluster.Execute();
    delete metric;
    delete refMatrix;
    const Rcpp::DataFrame clusterMetricsDataFrame = cluster.GetSensitivityData();
    const Rcpp::DataFrame iterationsMetricsDataFrame = cluster.GetClusterMetrics();
    const auto label = result->GetListVector().label;
    const Rcpp::DataFrame clusterDataFrame = result->GetListVector().listVector.CreateDataFrameFromList(
        featureColumnName, binColumnName);
    const Rcpp::DataFrame tidySharedDataFrame = CreateSharedDataFrame(combinedCountTableAdapter, result, binColumnName);
    delete(result);

    return Rcpp::List::create(Rcpp::Named("label") = std::stod(label),
      Rcpp::Named("abundance") = tidySharedDataFrame,
      Rcpp::Named("cluster") = clusterDataFrame,
      Rcpp::Named("cluster_metrics") = clusterMetricsDataFrame,
      Rcpp::Named("iteration_metrics") = iterationsMetricsDataFrame);
    return Rcpp::List::create();
}

//[[Rcpp::export]]
Rcpp::List OptiSplit(const SEXP& combinedData, const Rcpp::DataFrame& refList, const std::vector<std::string>& accnos,
    const float fitPercent, const std::string& featureColumnName, const std::string& binColumnName,
    const double cutoff, const bool isClosed = true, const bool printRef = false, const bool selfReference = true) {

    const ListVector refListOtuVector = Utils::CreateListVectorFromOtuList(refList["bin_name"],
          refList["sequence_name"]);
    const Rcpp::XPtr<DistanceFileReader> combinedDistanceData(combinedData);
    const CountTableAdapter combinedCountTableAdapter = combinedDistanceData.get()->GetCountTableAdapter();
    const SparseDistanceMatrix* combinedSparseMartix =  combinedDistanceData.get()->GetSparseMatrix();
    const ListVector* combinedListVector = combinedDistanceData.get()->GetListVector();
    const bool combinedIsSim = combinedDistanceData.get()->GetIsSimularity();
    const OptimatrixAdapter combinedOptiAdapter(cutoff);
    const auto* combinedOptiMatrix = combinedOptiAdapter.ConvertToOptimatrix(combinedSparseMartix, combinedListVector, combinedIsSim);
    delete combinedSparseMartix;
    delete combinedListVector;

    auto* splitMatrix = new SplitMatrix(combinedOptiMatrix, combinedCountTableAdapter,
        {accnos.begin(), accnos.end()});
    ClusterMetric* metric = new MCC();
    OptiFitCluster cluster(refMatrix, metric, refListOtuVector, cutoff, 0, selfReference, printRef, isClosed);
    const auto* result = cluster.Execute();
    delete metric;
    delete refMatrix;
    const Rcpp::DataFrame clusterMetricsDataFrame = cluster.GetSensitivityData();
    const Rcpp::DataFrame iterationsMetricsDataFrame = cluster.GetClusterMetrics();
    const auto label = result->GetListVector().label;
    const Rcpp::DataFrame clusterDataFrame = result->GetListVector().listVector.CreateDataFrameFromList(
        featureColumnName, binColumnName);
    const Rcpp::DataFrame tidySharedDataFrame = CreateSharedDataFrame(combinedCountTableAdapter, result, binColumnName);
    delete(result);
    return Rcpp::List::create(Rcpp::Named("label") = std::stod(label),
      Rcpp::Named("abundance") = tidySharedDataFrame,
      Rcpp::Named("cluster") = clusterDataFrame,
      Rcpp::Named("cluster_metrics") = clusterMetricsDataFrame,
      Rcpp::Named("iteration_metrics") = iterationsMetricsDataFrame);
    return Rcpp::List::create();
}
