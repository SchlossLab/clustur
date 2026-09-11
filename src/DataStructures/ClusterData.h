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

struct ClusterData {
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
        std::vector<ClusterExport*> results{};
        double cutoff{};

        ClusterData() = default;
        explicit ClusterData(const std::vector<OptiDataComponent>& dividedData,
            ClusterParameters parameters, ClusterMetric* metric,  const double cutoff) :
        dividedData(dividedData), clusterParameters(parameters), metric(metric), cutoff(cutoff){}
        // clusterData(bool showab, bool cla, bool df, vector< map<string, string> > dN, bool cns, double cu, int prec, int len, string meth, string opd, string vl, string ty) {
            // showabund = showab;
            // distNames = dN;
            // cutoff = cu;
            // classic = cla;
            // method = meth;
            // precision = prec;
            // length = len;
            // outputDir = opd;
            // vsearchLocation = vl;
            // deleteFiles = df;
            // cutoffNotSet = cns;
            // m = MothurOut::getInstance();
            // count = 0;
            // type = ty;
            // useName = false;
            // useCount = false;
        //     //numSingletons = 0;
        // }
        // void setOptiOptions(string metn, double stabMet, string init, int mxi ) {
        //     metricName = metn;
        //     stableMetric = stabMet;
        //     maxIters = mxi;
        //     initialize = init;
        // }
        // void setNamesCount(string cnf) {
        //     useName = false;
        //     useCount = false;
        //     if (type == "name") { useName = true;  }
        //     if (type == "count") { useCount = true; }
        // }
    };

#endif //REFACTOR_CLUSTERDATA_H