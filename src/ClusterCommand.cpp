/*
 *  clustercommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "MothurDependencies/ClusterCommand.h"
#include "MothurDependencies/AverageLinkage.h"
#include "MothurDependencies/Cluster.h"
#include "MothurDependencies/CompleteLinkage.h"
#include "MothurDependencies/OpticlusterData.h"
#include "MothurDependencies/SingleLinkage.h"
#include "MothurDependencies/WeightedLinkage.h"
#include "MothurDependencies/ClusterData.h"
#include "TestHelpers/TestHelper.h"
#include <string>

#include "Adapters/DataFrameAdapter.h"
using namespace std;

ClusterCommand::~ClusterCommand() {
}

// //**********************************************************************************************************************
/// Bad allocations, returns basic_string, returns empty string, returns non-utf8 characters, etc
/// @param optiMatrix
/// @param cutoffValue
/// @return
ClusterExport* ClusterCommand::runOptiCluster(OptiMatrix *optiMatrix, const double cutoffValue) {
    cutoff = cutoffValue;
    // clusterMetrics += ("\nClustering " + distfile + "\n");
    initialize = "singleton";
    auto* data = new OpticlusterData("", cutoff);
    ClusterMetric *metric = nullptr;
    std::string clusterMetrics;
    std::string sensFile;
    metricName = "mcc";
    if (metricName == "mcc") { metric = new MCC(); } else if (
        metricName == "sens") { metric = new Sensitivity(); } else if (metricName == "spec") {
        metric = new Specificity();
    } else if (
        metricName == "tptn") { metric = new TPTN(); } else if (metricName == "tp") { metric = new TP(); } else if (
        metricName == "tn") { metric = new TN(); } else if (metricName == "fp") { metric = new FP(); } else if (
        metricName == "fn") { metric = new FN(); } else if (
        metricName == "f1score") { metric = new F1Score(); } else if (
        metricName == "accuracy") { metric = new Accuracy(); } else if (
        metricName == "ppv") { metric = new PPV(); } else if (metricName == "npv") { metric = new NPV(); } else if (
        metricName == "fdr") { metric = new FDR(); } else if (metricName == "fpfn") { metric = new FPFN(); } else {
        return {};
    }
    // Setting headers


    std::vector<std::string> sensfileHeaders{"label","cutoff","ttp","tn","fp","fn","sensitivity",
        "specificity","ppv","npv","fdr","accuracy","mcc","f1score"};
    // sensFile = "label\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n";

    std::vector<std::string> clusterMetricsHeaders{"iter","time","label","num_otus","cutoff","tp","tn",
      "fp","fn","sensitivity","specificity","ppv","npv",
        "fdr", "accuracy", "mcc", "f1score"};
    // clusterMetrics = (
    //     "iter\ttime\tlabel\tnum_otus\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n");
    bool printHeaders = true;
    cutoffs.insert(std::to_string(cutoff));

    OptiData *matrix = optiMatrix;
    for (auto it = cutoffs.begin(); it != cutoffs.end(); it++) {
        OptiCluster cluster(matrix, metric, 0);
        int iters = 0;
        double listVectorMetric = 0; //worst state
        double delta = 1;
        long long numBins;
        double tp, tn, fp, fn;
        vector<double> stats;
        std::vector<std::string> clusterMetricList;
        cluster.initialize(listVectorMetric, canShuffle, initialize);
        stats = cluster.getStats(tp, tn, fp, fn);
        numBins = cluster.getNumBins();
        clusterMetrics = ("0,0," + std::to_string(cutoff) + "," + std::to_string(numBins) + "," +
                           std::to_string(cutoff) + "," + std::to_string(tp) + "," + std::to_string(tn) + "," +
                           std::to_string(fp) + "," + std::to_string(fn) + ",");


        for (double result: stats) {
            clusterMetrics += (std::to_string(result) + ",");
        }
        util.AddRowToDataFrameMap(dataframeMapClusterMetrics, clusterMetrics, clusterMetricsHeaders);
        //m->mothurOutEndLine();
        // Stable Metric -> Keep the data stable, to prevent errors (rounding errors)
        // The difference between what the current and last metric (delta)
        // MaxIters -> is an exit condition
        while ((delta > stableMetric) && (iters < maxIters)) {
            //long start = std::time(nullptr);
            double oldMetric = listVectorMetric;

            cluster.update(listVectorMetric);

            delta = std::abs(oldMetric - listVectorMetric);
            iters++;

            stats = cluster.getStats(tp, tn, fp, fn);

            numBins = cluster.getNumBins();

            clusterMetrics = (std::to_string(iters) + "," + std::to_string(std::time(nullptr) - start) + "," +
                               std::to_string(cutoff) + "," + std::to_string(numBins) + "," +
                               std::to_string(cutoff) + "," + std::to_string(tp) + "," + std::to_string(tn) + ","
                               + std::to_string(fp) + "," + std::to_string(fn) + ",");
            for (double result: stats) {
                clusterMetrics += (std::to_string(result) + ",");
            }
            util.AddRowToDataFrameMap(dataframeMapClusterMetrics, clusterMetrics, clusterMetricsHeaders);
        }
        ListVector *list = nullptr;
        // clusterMetrics += "\n\n";
        list = cluster.getList();
        list->setLabel(std::to_string(cutoff));
        //
        if (printHeaders) {
            //only print headers the first time
            printHeaders = false;
        } else { list->setPrintedLabels(printHeaders); }
        OptiClusterInformation clusterInformation;
        clusterInformation.label = std::to_string(cutoff);
        clusterInformation.numberOfOtu = static_cast<int>(numBins);
        clusterInformation.clusterBins = list->print(listFile);
        data->AddToData(clusterInformation);
        data->SetListVector(*list, std::to_string(cutoff));
        // clusterMatrixOutput = list->print(listFile);
        stats = cluster.getStats(tp, tn, fp, fn);

        sensFile += std::to_string(cutoff) + ',' + std::to_string(cutoff) + ',' + std::to_string(tp) + ',' +
                std::to_string(tn) + ',' +
                std::to_string(fp) + ',' + std::to_string(fn) + ',';
        for (double result: stats) { sensFile += std::to_string(result) + ','; }
        util.AddRowToDataFrameMap(dataframeMapSensMetrics, sensFile, sensfileHeaders);
    }
    delete matrix;
    return data;
}

ClusterExport* ClusterCommand::runMothurCluster(const std::string &clusterMethod, SparseDistanceMatrix *matrix,
                                             double cutoff, ListVector *list) {
    //
    Cluster *cluster = nullptr;
    method = clusterMethod;
    // auto* result = new ClusterResult();
    auto* clusterData = new ClusterData("");
    auto rAbund = list->getRAbundVector();
    if (clusterMethod ==  "furthest")	{	cluster = new CompleteLinkage(&rAbund, list, matrix, cutoff, method, adjust); }
    else if(clusterMethod == "nearest"){	cluster = new SingleLinkage(&rAbund, list, matrix, cutoff, method, adjust); }
    else if(clusterMethod == "average"){	cluster = new AverageLinkage(&rAbund, list, matrix, cutoff, method, adjust);	}
    else {	cluster = new WeightedLinkage(&rAbund, list, matrix, cutoff, method, adjust);}

    map<string, int> counts;
    this->cutoff = cutoff;
    float previousDist = 0.00000;
    float rndPreviousDist = 0.00000;
    oldList = *list;
    bool printHeaders = false;
    std::string clusterResult;
    double highestDistLabel =  -1;
    std::string binResults;
    list->setPrintedLabels(printHeaders);
    while ((matrix->getSmallDist() <= cutoff) && (matrix->getNNodes() > 0)) { //TODO We are getting values that are just barely grater than 0, we need to figure out how to deal with them
        cluster->update(cutoff);
        ClusterInformation data;
        const float dist = matrix->getSmallDist(); // Round to the third decimal place
        //Rcpp::Rcout << dist << std::endl;
        const float rndDist = util.ceilDist(dist, precision);
        if (previousDist <= 0.0000 && !util.isEqual(dist, previousDist)) {
            data.label = "0.00000";
            data.numberOfOtu = oldList.getNumBins();
        } else if (!util.isEqual(rndDist, rndPreviousDist)) {
            data.label = std::to_string(rndPreviousDist);
            data.numberOfOtu = oldList.getNumBins();
        }
       
        
       
        if(!data.label.empty()) {
            data.clusterBins = oldList.print(listFile);
            auto* vec = new ListVector(oldList);
            list->setPrintedLabels(false);
            clusterData->AddToData(data);
            if(rndPreviousDist > highestDistLabel) {
                highestDistLabel = rndPreviousDist;
                vec->setLabel(std::to_string(highestDistLabel));
                clusterData->SetListVector(*vec, std::to_string(highestDistLabel)); // vec might be a shallow copy
            }
        }
        oldList = *list;
        previousDist = dist;
        rndPreviousDist = rndDist;

    }
    ClusterInformation data;
    if(previousDist <= 0.0000) {
        data.label = std::to_string(previousDist);
        data.numberOfOtu = oldList.getNumBins();
    }
    else if(rndPreviousDist<cutoff) {
        data.label = std::to_string(rndPreviousDist);
        data.numberOfOtu = oldList.getNumBins();
    }

    if(!data.label.empty()) {
        data.clusterBins = oldList.print(listFile);
        auto* vec = new ListVector(oldList);
        
        clusterData->AddToData(data);
        if(rndPreviousDist > highestDistLabel) {
            highestDistLabel = rndPreviousDist;
            vec->setLabel(std::to_string(highestDistLabel));
            clusterData->SetListVector(*vec, std::to_string(highestDistLabel));
        }
    }
    delete(cluster);
    return clusterData;
}

std::string ClusterCommand::PrintData(const string& label, map<string, int> &counts, bool &ph) {
    //TODO Return a ClusterData to add to cluster results here!
    oldList.setPrintedLabels(ph);
    ph = false;
    oldList.setLabel(label);
    std::string data = label + "\t" + std::to_string(oldList.getNumBins());
    if (!countfile.empty()) {
        data += oldList.print(listFile, counts);
    } else {
        data += oldList.print(listFile);
    }
    return data;
}


