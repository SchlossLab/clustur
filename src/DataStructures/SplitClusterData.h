//
// Created by Gregory Johnson on 9/9/26.
//

#ifndef REFACTOR_CLUSTERDATA_H
#define REFACTOR_CLUSTERDATA_H
#include "../Clusters/ClusterMethod.h"
#include "../DataStructures/OptiData.h"
#include <vector>

#include "ClusterParameters.h"
#include "OptidataComponents.h"

struct SplitClusterData {
    // Utils util;
    // int count, precision, length, maxIters; //numSingletons,
    // bool showabund, classic, useName, useCount, deleteFiles, cutoffNotSet;
    // double cutoff, stableMetric;
    // ofstream outList, outRabund, outSabund;
    // string tag, method,  vsearchLocation, metricName, initialize, outputDir, type;
    // vector< map<string, string> > distNames;
    // set<string> labels;
    // vector<string> listFileNames;
    std::vector<OptiDataComponent> dividedData;
    ClusterParameters clusterParameters{};
    ClusterMetric* metric;
    RandomNumberSitmo rng;
    std::vector<ClusterExport*> results{};
    double cutoff{};

    SplitClusterData() = default;
    explicit SplitClusterData(const std::vector<OptiDataComponent>& dividedData,
        ClusterParameters parameters, ClusterMetric* metric, const RandomNumberSitmo& rng, const double cutoff) :
    dividedData(dividedData), clusterParameters(parameters), metric(metric), rng(rng), cutoff(cutoff){}
    // ~ClusterData() {
    //     for (auto& res : results) {
    //         delete res;
    //     }
    // }

};

#endif //REFACTOR_CLUSTERDATA_H