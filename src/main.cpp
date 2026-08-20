#include <Rcpp.h>
#include <cctype>
#include <algorithm>

#include "Adapters/CountTableAdapter.h"
#include "Adapters/MatrixAdapter.h"
#include "Adapters/OptimatrixAdapter.h"
#include "Clusters/AverageLinkage.h"
#include "Clusters/CompleteLinkage.h"
#include "Clusters/SingleLinkage.h"
#include "Clusters/WeightedLinkage.h"
#include "Clusters/OptiCluster.h"
#include "Clusters/Optifitcluster.h"
#include "Clusters/Metrics/mcc.h"
#include "Clusters/Metrics/tptn.h"
#include "DataExporters/ClusterExport.h"
#include "SharedFileData/SharedFileBuilder.h"
#include "FileReaders/DistanceFileReader.h"
#include "FileReaders/ColumnDistanceMatrixReader.h"
#include "FileReaders/ReadPhylipMatrix.h"
#include "MothurDependencies/CreateDistanceMatrix.h"
#include "MothurDependencies/OneGapPairwiseDistance.h"
#include "MothurDependencies/PairwiseDistanceCalculator.h"


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
    CountTableAdapter countTableAdapter;
    countTableAdapter.CreateDataFrameMap(countTable);
    if(isPhylip) {
        DistanceFileReader* read = new ReadPhylipMatrix(cutoff, isSim);
        read->SetCountTableAdapter(countTableAdapter);
        if (const bool succeeded = read->Read(filePath); !succeeded) {
            const std::set<std::string> params = read->GetFailureParameters();
            delete read;

            Utils::CheckForDistanceFileError(params);
        }

        return Rcpp::XPtr<DistanceFileReader>(read);
    }
    DistanceFileReader* read = new ColumnDistanceMatrixReader(cutoff, isSim);
    read->SetCountTableAdapter(countTableAdapter);
    if (const bool succeeded = read->Read(filePath); !succeeded) {
        const std::set<std::string> params = read->GetFailureParameters();
        delete read;
        Utils::CheckForDistanceFileError(params);
    }
    return Rcpp::XPtr<DistanceFileReader>(read);
}

 //[[Rcpp::export]]
 SEXP ProcessSparseMatrix(const std::vector<int> &xPosition,
     const std::vector<int> &yPosition, const std::vector<double> &data, const Rcpp::DataFrame& countTable,
     const double cutoff, const bool isSim) {
     CountTableAdapter countTableAdapter;
     countTableAdapter.CreateDataFrameMap(countTable);
     const MatrixAdapter adapter(xPosition, yPosition, data, cutoff, isSim, countTableAdapter);
     auto* read = new DistanceFileReader(adapter.CreateSparseMatrix(),
         adapter.CreateListVector(), countTableAdapter, cutoff, isSim);
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

ClusterMethod* GetClusterMethod(const std::string& method, ListVector* listVector,
    SparseDistanceMatrix* matrix, RAbundVector& rAbund, const double cutoff, const double adjust = -1) {
    if (method == "furthest")	return new CompleteLinkage(&rAbund, listVector, matrix, cutoff, method, adjust);
    if(method == "nearest") return new SingleLinkage(&rAbund, listVector, matrix, cutoff, method, adjust);
    if(method == "average")	return new AverageLinkage(&rAbund, listVector, matrix, cutoff, method, adjust);
    return new WeightedLinkage(&rAbund, listVector, matrix, cutoff, method, adjust);
}


ListVector CreateListVectorFromOtuList(const std::vector<std::string> &otuBins, const std::vector<std::string> &sequences) {
    const int size = otuBins.size();
    std::string lastOtu = otuBins[otuBins.size() - 1];
    lastOtu.erase(std::remove_if(lastOtu.begin(), lastOtu.end(), isalpha), lastOtu.end());
    const int otuSize = std::stoi(lastOtu) - 1;
    ListVector listVector(otuSize);
    std::string currentOtuBin = otuBins[0];
    std::string currentBinData = "";
    size_t index = 0;
    for (int i = 0; i < size; i++) {
        if (currentOtuBin == otuBins[i]) {
            currentBinData.append(sequences[i] + ",");
            continue;
        }
        currentOtuBin = otuBins[i];
        currentBinData.append(sequences[i] + ",");
        listVector.set(index++, {currentBinData.begin(), currentBinData.end() - 1});
        currentBinData = "";
    }
    return listVector;
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
    ClusterMethod* clusterMethod = GetClusterMethod(method, listVector, sparseMatrix,
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
    const double cutoff, const double fitPercent = 50) {
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
    OptiFitCluster cluster(refMatrix, metric,"denovo", cutoff, 0);
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
    const double cutoff) {
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
    OptiFitCluster cluster(refMatrix, metric,"userref", cutoff, 0);

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
Rcpp::List OptiFit3(const SEXP& combinedData, const Rcpp::DataFrame& refList, const float fitPercent, const std::string& featureColumnName, const std::string& binColumnName,
    const double cutoff) {
// fitPercent = fitPercent = ((count-refCount) / static_cast<float>(count));
    const ListVector refListOtuVector = CreateListVectorFromOtuList(refList["bin_name"],
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
    //
    // const Rcpp::XPtr<DistanceFileReader> fitDistanceData(fitData);
    // const CountTableAdapter fitCountTableAdapter = fitDistanceData.get()->GetCountTableAdapter();
    // const SparseDistanceMatrix* fitSparseMatrix =  fitDistanceData.get()->GetSparseMatrix();
    // const ListVector* fitListVector = fitDistanceData.get()->GetListVector();
    // const bool fitIsSim = fitDistanceData.get()->GetIsSimularity();
    // const OptimatrixAdapter fitOptiAdapter(cutoff);
    // const auto* fitOptiMatrix = fitOptiAdapter.ConvertToOptimatrix(fitSparseMatrix, fitListVector, fitIsSim);
    // FastaDatabase fitFastaDatabase(fitFasta["sequence_name"], fitFasta["sequence"]);
    // delete fitSparseMatrix;
    // delete fitListVector;




    auto* refMatrix = new OptiRefMatrix(combinedOptiMatrix, combinedCountTableAdapter, fitPercent, cutoff);
    ClusterMetric* metric = new MCC();
    OptiFitCluster cluster(refMatrix, metric, refListOtuVector, cutoff, 0);
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
void ToListVector(const Rcpp::DataFrame& df) {
   ListVector vec = CreateListVectorFromOtuList(df["bin_name"], df["sequence_name"]);
   size_t v = vec.size();
}

//[[Rcpp::export]]
Rcpp::DataFrame CreateDataFrameFromSparseCountTable(const Rcpp::DataFrame& countTable) {
    CountTableAdapter adapter;
    adapter.CreateDataFrameMapFromSparseCountTable(countTable);
    return adapter.ReCreateDataFrame();
}

//[[Rcpp::export]]
double GetDist(const std::string& sequenceOne, const std::string& sequenceTwo) {
    PairwiseDistanceCalculator* calculator = new OneGapPairwiseDistance();
    return calculator->Execute(sequenceOne, sequenceTwo);
}


//[[Rcpp::export]]
double CreateSparseMatrix(const std::vector<std::string>& sequences, const double cutoff) {
    PairwiseDistanceCalculator* calculator = new OneGapPairwiseDistance();
    SparseDistanceMatrix matrix = CreateDistanceMatrix::CreateSparseDistanceMatrix(sequences, calculator, cutoff);
    delete calculator;
    return 0;

}


//[[Rcpp::export]]
double AddDataToDistanceData(SEXP& refData, const SEXP& fitData,
    const Rcpp::DataFrame& refFasta,
    const Rcpp::DataFrame& fitFasta,
    const double cutoff) {

    const Rcpp::XPtr<DistanceFileReader> refDistanceData(refData);
    const CountTableAdapter refCountTableAdapter = refDistanceData.get()->GetCountTableAdapter();
    const FastaDatabase refFastaDatabase(refFasta["sequence_name"], refFasta["sequence"]);

    const Rcpp::XPtr<DistanceFileReader> fitDistanceData(fitData);
    const CountTableAdapter fitCountTableAdapter = fitDistanceData.get()->GetCountTableAdapter();
    const SparseDistanceMatrix* fitSparseMatrix =  fitDistanceData.get()->GetSparseMatrix();
    const ListVector* fitListVector = fitDistanceData.get()->GetListVector();
    FastaDatabase fitFastaDatabase(fitFasta["sequence_name"], fitFasta["sequence"]);

    refDistanceData.get()->AddFittedDataToReference(fitSparseMatrix, fitListVector,
        fitCountTableAdapter, refFastaDatabase, fitFastaDatabase, cutoff);
    const double fitCountSize = static_cast<double>(fitCountTableAdapter.GetSequences().size());
    const double refCountSize = static_cast<double>(refCountTableAdapter.GetSequences().size());
    const double result = fitCountSize / (fitCountSize + refCountSize);
    return result;// Fit Percentage
}

//[[Rcpp::export]]
SEXP CopyObject(const SEXP& distanceObject) {
    const Rcpp::XPtr<DistanceFileReader> refDistanceData(distanceObject);
    DistanceFileReader* copy = new DistanceFileReader(*refDistanceData.get());
    return Rcpp::XPtr<DistanceFileReader>(copy);
}