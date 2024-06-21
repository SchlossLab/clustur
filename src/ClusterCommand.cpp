/*
 *  clustercommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "MothurDependencies/ClusterCommand.h"

#include "MothurDependencies/Cluster.h"
#include "MothurDependencies/SingleLinkage.h"
using namespace std;
ClusterCommand::~ClusterCommand() {
}

// //**********************************************************************************************************************
/// Bad allocations, returns basic_string, returns empty string, returns non-utf8 characters, etc
/// @param optiMatrix
/// @return
std::vector<std::string> ClusterCommand::runOptiCluster(OptiMatrix *optiMatrix) {
    std::string clusterMetrics;
    std::string sensFile;
    std::string outStep;
    std::string clusterMatrixOutput;
    if (!cutOffSet) {
        // clusterMetrics += ("\nYou did not set a cutoff, using 0.03.\n");
        cutoff = 0.05;
    }

    // clusterMetrics += ("\nClustering " + distfile + "\n");

    ClusterMetric *metric = nullptr;
    metricName = "mcc";
    if (metricName == "mcc") { metric = new MCC(); } else if (metricName == "sens") { metric = new Sensitivity(); } else
        if (metricName == "spec") { metric = new Specificity(); } else if (
            metricName == "tptn") { metric = new TPTN(); } else if (metricName == "tp") { metric = new TP(); } else if (
            metricName == "tn") { metric = new TN(); } else if (metricName == "fp") { metric = new FP(); } else if (
            metricName == "fn") { metric = new FN(); } else if (
            metricName == "f1score") { metric = new F1Score(); } else if (
            metricName == "accuracy") { metric = new Accuracy(); } else if (
            metricName == "ppv") { metric = new PPV(); } else if (metricName == "npv") { metric = new NPV(); } else if (
            metricName == "fdr") { metric = new FDR(); } else if (metricName == "fpfn") { metric = new FPFN(); } else {
            return {};
        }

    // string nameOrCount = "";
    // string thisNamefile = "";
    // map<string, int> counts;
    // if (countfile != "") {
    //     nameOrCount = "count";
    //     thisNamefile = countfile;
    // } //CountTable ct; ct.readTable(countfile, false, false); counts = ct.getNameMap(); }
    // else if (namefile != "") {
    //     nameOrCount = "name";
    //     thisNamefile = namefile;
    // }

    // string distfile = columnfile;
    // if (format == "phylip") { distfile = phylipfile; }
    //
    sensFile += "label\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n";
    clusterMetrics += (
        "iter\ttime\tlabel\tnum_otus\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n");
    outStep +=
            "iter\ttime\tlabel\tnum_otus\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n";

    bool printHeaders = true;

    cutoffs.insert("0.2");

    OptiData *matrix = optiMatrix;
    for (auto it = cutoffs.begin(); it != cutoffs.end(); it++) {
        OptiCluster cluster(matrix, metric, 0);

        int iters = 0;
        double listVectorMetric = 0; //worst state
        double delta = 1;
        long long numBins;
        double tp, tn, fp, fn;
        vector<double> stats;
        cluster.initialize(listVectorMetric, canShuffle, initialize);
        stats = cluster.getStats(tp, tn, fp, fn);
        numBins = cluster.getNumBins();
        clusterMetrics += ("0\t0\t" + std::to_string(cutoff) + "\t" + std::to_string(numBins) + "\t" +
                              std::to_string(cutoff) + "\t" + std::to_string(tp) + "\t" + std::to_string(tn) + "\t" +
                              std::to_string(fp) + "\t" + std::to_string(fn) + "\t");
        outStep += "0\t0\t" + std::to_string(cutoff) + "\t" + std::to_string(numBins) + "\t" + std::to_string(cutoff) +
                "\t" + std::to_string(tp) + '\t' + std::to_string(tn) + '\t' + std::to_string(fp) + '\t' +
                std::to_string(fn) + '\t';
        for (double result : stats) {
            clusterMetrics += (std::to_string(result) + "\t");
            outStep += std::to_string(result) + "\t";
        }
        //m->mothurOutEndLine();
        clusterMetrics += "\n";
        outStep += "\n";
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

            clusterMetrics += (std::to_string(iters) + "\t" + std::to_string(time(nullptr) - start) + "\t" +
                                  std::to_string(cutoff) + "\t" + std::to_string(numBins) + "\t" +
                                  std::to_string(cutoff) + "\t" + std::to_string(tp) + "\t" + std::to_string(tn) + "\t"
                                  + std::to_string(fp) + "\t" + std::to_string(fn) + "\t");
            outStep += (std::to_string(iters) + "\t" + std::to_string(time(nullptr) - start) + "\t" +
                        std::to_string(cutoff) + "\t" + std::to_string(numBins) + "\t" + std::to_string(cutoff) + "\t")
                    + std::to_string(tp) + '\t' + std::to_string(tn) + '\t' + std::to_string(fp) + '\t' +
                    std::to_string(fn) +
                    '\t';
            for (double result : stats) {
                clusterMetrics += (std::to_string(result) + "\t");
                outStep += std::to_string(result) + "\t";
            }
            clusterMetrics += "\n";
            outStep += "\n";
        }
        ListVector *list = nullptr;

        // clusterMetrics += "\n\n";
        list = cluster.getList();
        //
        if (printHeaders) {
            //only print headers the first time
            printHeaders = false;
        } else { list->setPrintedLabels(printHeaders); }
        clusterMatrixOutput = list->print(listFile);
        delete list;
        stats = cluster.getStats(tp, tn, fp, fn);

        sensFile += std::to_string(cutoff) + '\t' + std::to_string(cutoff) + '\t' + std::to_string(tp) + '\t' +
                std::to_string(tn) + '\t' +
                std::to_string(fp) + '\t' + std::to_string(fn) + '\t';
        for (double result : stats) { sensFile += std::to_string(result) + '\t'; }
     }
    delete matrix;
    return {clusterMatrixOutput, sensFile, clusterMetrics};
}

int ClusterCommand::runMothurCluster(SingleLinkage* cluster, SparseDistanceMatrix* matrix, double cutoff, const ListVector* list){
        //
        map<string, int> counts;
        this->cutoff = cutoff;
        float previousDist = 0.00000;
        float rndPreviousDist = 0.00000;
        oldList = *list;
        print_start = true;
        start = time(nullptr);
        loops = 0;
        double saveCutoff = cutoff;
        bool printHeaders = true;

    while ((matrix->getSmallDist() <= cutoff) && (matrix->getNNodes() > 0)){

        // if (m->getControl_pressed()) { //clean up
        //     delete list; delete matrix; delete rabund; delete cluster;
        //     if(countfile == "") {rabundFile.close(); sabundFile.close();  util.mothurRemove((fileroot+ tag + ".rabund")); util.mothurRemove((fileroot+ tag + ".sabund")); }
        //     listFile.close(); util.mothurRemove((fileroot+ tag + ".list")); outputTypes.clear(); return 0;
        // }
        //
        // if (print_start && util.isTrue(timing)) {
        //     m->mothurOut("Clustering (" + tag + ") dist " + toString(matrix->getSmallDist()) + "/"
        //                  + toString(util.roundDist(matrix->getSmallDist(), precision))
        //                  + "\t(precision: " + toString(precision) + ", Nodes: " + toString(matrix->getNNodes()) + ")");
        //     cout.flush();
        //     print_start = false;
        // }

        cluster->update(cutoff);
        std::cout << cutoff << std::endl;
        const double smallest_dist =  matrix->getSmallDist();
        const int nNodes = matrix -> getNNodes();

        //std::cout << "Nodes: " << nNodes << "\t" << "smallest_dist: " << smallest_dist << std::endl;

        float dist = matrix->getSmallDist();
        float rndDist = util.ceilDist(dist, precision);

        if(previousDist <= 0.0000 && !util.isEqual(dist, previousDist)) {
            printData("unique", counts, printHeaders);
        }
        else if(!util.isEqual(rndDist, rndPreviousDist)) {
            printData(std::to_string(rndPreviousDist), counts, printHeaders);
        }
        oldList = *list;
        previousDist = dist;
        rndPreviousDist = rndDist;
    }
     //const std::string data = oldList.print(listFile);
     //std::cout << "\n" <<data;



        return 0;
}
void ClusterCommand::printData(const string label, map<string, int>& counts, bool& ph){
        oldList.setPrintedLabels(ph); ph = false;
        print_start = true;
        loops = 0;
        start = time(nullptr);

        oldList.setLabel(label);
    std::string data;
        if(countfile != "") {
            data = oldList.print(listFile, counts);
        }else {
            data = oldList.print(listFile);
        }
    std::cout << data;
}