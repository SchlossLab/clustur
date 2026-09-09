//
// Created by Gregory Johnson on 9/2/26.
//

#ifndef REFACTOR_CLUSTERSPLIT_H
#define REFACTOR_CLUSTERSPLIT_H
#include "ClusterMethod.h"
#include "../DataStructures/FastaDatabase.h"
#include "../DataStructures/TaxonomyData.h"
#include "../MothurDependencies/PairwiseDistanceCalculator.h"


class ClusterSplit : ClusterMethod {
public:
    ClusterSplit(const FastaDatabase& fastaDatabase, const std::vector<TaxonomyData>& taxaData,
        PairwiseDistanceCalculator* calculator, double cutoff, int taxonomyCutoff);
    ~ClusterSplit() override;
    ClusterExport * Execute() override;
private:
    FastaDatabase fastaData;
    PairwiseDistanceCalculator* calculator;
    std::vector<TaxonomyData> taxaData;
    int taxonomyCutoff;
    double cutoff;

    struct clusterData {
        Utils util;
        int count, precision, length, maxIters; //numSingletons,
        bool showabund, classic, useName, useCount, deleteFiles, cutoffNotSet;
        double cutoff, stableMetric;
        ofstream outList, outRabund, outSabund;
        string tag, method,  vsearchLocation, metricName, initialize, outputDir, type;
        vector< map<string, string> > distNames;
        set<string> labels;
        vector<string> listFileNames;

        clusterData(){}
        clusterData(bool showab, bool cla, bool df, vector< map<string, string> > dN, bool cns, double cu, int prec, int len, string meth, string opd, string vl, string ty) {
            showabund = showab;
            distNames = dN;
            cutoff = cu;
            classic = cla;
            method = meth;
            precision = prec;
            length = len;
            outputDir = opd;
            vsearchLocation = vl;
            deleteFiles = df;
            cutoffNotSet = cns;
            m = MothurOut::getInstance();
            count = 0;
            type = ty;
            useName = false;
            useCount = false;
            //numSingletons = 0;
        }
        void setOptiOptions(string metn, double stabMet, string init, int mxi ) {
            metricName = metn;
            stableMetric = stabMet;
            maxIters = mxi;
            initialize = init;
        }
        void setNamesCount(string cnf) {
            useName = false;
            useCount = false;
            if (type == "name") { useName = true;  }
            if (type == "count") { useCount = true; }
        }
    };
};


#endif //REFACTOR_CLUSTERSPLIT_H