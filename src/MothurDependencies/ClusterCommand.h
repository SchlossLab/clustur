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
	std::string runOptiCluster(OptiMatrix*);


	
private:
	ListVector* list;
	ListVector oldList;
	bool abort, sim, cutOffSet;
	string method, fileroot, tag, phylipfile, columnfile, namefile, format, distfile, countfile, fastafile, inputDir, vsearchLocation, metric, initialize;
	double cutoff, stableMetric = 0;
    float adjust;
	string showabund, timing, metricName;
	int precision, length, maxIters = 2, processors;
	ofstream sabundFile, rabundFile, listFile;
    set<string> cutoffs;
	Utils util;
	bool print_start;
	time_t start;
	unsigned long loops;
	vector<string> outputNames;
};

#endif
