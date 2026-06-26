#ifndef CLUSTERCOMMAND_H
#define CLUSTERCOMMAND_H
/*
 *  clustercommand.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include <vector>
#include <string>

#include "OptiCluster.h"
#include "OptiData.h"
#include "Utils.h"
#include <set>
#include <fstream>
#include <map>
#include "OptiMatrix.h"
#include "ListVector.h"
#include <chrono>
#include "Metrics/accuracy.hpp"
#include "Metrics/f1score.hpp"
#include "Metrics/fdr.hpp"
#include "Metrics/fn.hpp"
#include "Metrics/fp.hpp"
#include "Metrics/fpfn.hpp"
#include "Metrics/mcc.hpp"
#include "Metrics/npv.hpp"
#include "Metrics/ppv.hpp"
#include "Metrics/sensitivity.hpp"
#include "Metrics/specificity.hpp"
#include "Metrics/tn.hpp"
#include "Metrics/tp.hpp"
#include "Metrics/tptn.hpp"
#include "ClusterMetric.h"
#include "OptiMatrix.h"
#include <chrono>
#include <Rcpp.h>

#include "ClusterExport.h"
#include "SingleLinkage.h"
#include "../Adapters/DataFrameAdapter.h"
using namespace std;

/* The cluster() command:
	The cluster command outputs a .list , .rabund and .sabund files.  
	The cluster command parameter options are method, cuttoff and precision. No parameters are required.  
	The cluster command should be in the following format: cluster(method=yourMethod, cutoff=yourCutoff, precision=yourPrecision).  
	The acceptable methods are furthest, nearest and average.  If you do not provide a method the default algorithm is furthest neighbor.  
	The cluster() command outputs three files *.list, *.rabund, and *.sabund.   */


class ClusterCommand {
public:
	//ClusterCommand(string);
	ClusterCommand() {}
	~ClusterCommand();
	bool SetMaxIterations(const int iterations) {maxIters = iterations; return maxIters == iterations;}
	bool SetOpticlusterRandomShuffle(const bool shuffle) {canShuffle = shuffle; return canShuffle;}
	bool SetMetricType(const string& newMetric) {metric = newMetric; return metric == newMetric;}
	Rcpp::DataFrame GetSensitivityData() const {return DataFrameAdapter::UnorderedMapToDataFrame(dataframeMapSensMetrics);}
	Rcpp::DataFrame GetClusterMetrics() const {return DataFrameAdapter::UnorderedMapToDataFrame(dataframeMapClusterMetrics);}
	ClusterExport* runOptiCluster(OptiMatrix*, double);

	ClusterExport* runMothurCluster(const std::string& clusterMethod, SparseDistanceMatrix *matrix, double cutoff, ListVector*);

	std::string PrintData(const string &label, map<string, int> &counts, bool &ph);


private:
	ListVector oldList;
	string method, fileroot, tag, phylipfile, columnfile, namefile, format, distfile, countfile, fastafile, inputDir, vsearchLocation, metric, initialize;
	double cutoff, stableMetric = 0;
    float adjust = -1;
	string showabund, timing, metricName;
	int precision = 100, maxIters = 100;
	ofstream sabundFile, rabundFile, listFile;
    set<string> cutoffs;
	Utils util;
	bool canShuffle = true;
	vector<string> outputNames;
	std::unordered_map<std::string, std::vector<std::string>> dataframeMapClusterMetrics;
	std::unordered_map<std::string, std::vector<std::string>> dataframeMapSensMetrics;
};

#endif
