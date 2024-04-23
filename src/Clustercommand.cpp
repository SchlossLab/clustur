/*
 *  clustercommand.cpp
 *  Dotur
 *
 *  Created by Sarah Westcott on 1/2/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "MothurDependencies/Clustercommand.h"
using namespace std;
ClusterCommand::~ClusterCommand() {
}

// //**********************************************************************************************************************
/// Bad alloctions, returns basic_string, returns emptry string, returns non-utf8 characters, etc
/// @param optiMatrix
/// @return
int ClusterCommand::runOptiCluster(OptiMatrix *optiMatrix) {
    if (!cutOffSet) {
        clusterOutputData += ("\nYou did not set a cutoff, using 0.03.\n");
        cutoff = 0.5;
    }

    clusterOutputData += ("\nClustering " + distfile + "\n");

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
            return 0;
        }

    string nameOrCount = "";
    string thisNamefile = "";
    map<string, int> counts;
    if (countfile != "") {
        nameOrCount = "count";
        thisNamefile = countfile;
    } //CountTable ct; ct.readTable(countfile, false, false); counts = ct.getNameMap(); }
    else if (namefile != "") {
        nameOrCount = "name";
        thisNamefile = namefile;
    }

    string distfile = columnfile;
    if (format == "phylip") { distfile = phylipfile; }

    sensFile += "label\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n";
    clusterOutputData += (
        "\n\niter\ttime\tlabel\tnum_otus\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n");
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
        vector<double> results;
        cluster.initialize(listVectorMetric, true, initialize);
        results = cluster.getStats(tp, tn, fp, fn);
        numBins = cluster.getNumBins();
        clusterOutputData += ("0\t0\t" + std::to_string(cutoff) + "\t" + std::to_string(numBins) + "\t" +
                              std::to_string(cutoff) + "\t" + std::to_string(tp) + "\t" + std::to_string(tn) + "\t" +
                              std::to_string(fp) + "\t" + std::to_string(fn) + "\t");
        outStep += "0\t0\t" + std::to_string(cutoff) + "\t" + std::to_string(numBins) + "\t" + std::to_string(cutoff) +
                "\t" + std::to_string(tp) + '\t' + std::to_string(tn) + '\t' + std::to_string(fp) + '\t' +
                std::to_string(fn) + '\t';
        for (double result : results) {
            clusterOutputData += (std::to_string(result) + "\t");
            outStep += std::to_string(result) + "\t";
        }
        //m->mothurOutEndLine();
        outStep += "\n";
        // Stable Metric -> Keep the data stable, to prevent errors (rounding errors)
        // The difference between what the current and last metric (delta)
        // MaxIters -> is an exit condition
        while ((delta > stableMetric) && (iters < maxIters)) {
            long start = time(nullptr);
            double oldMetric = listVectorMetric;

            cluster.update(listVectorMetric);

            delta = std::abs(oldMetric - listVectorMetric);
            iters++;

            results = cluster.getStats(tp, tn, fp, fn);
            numBins = cluster.getNumBins();

            clusterOutputData += (std::to_string(iters) + "\t" + std::to_string(time(nullptr) - start) + "\t" +
                                  std::to_string(cutoff) + "\t" + std::to_string(numBins) + "\t" +
                                  std::to_string(cutoff) + "\t" + std::to_string(tp) + "\t" + std::to_string(tn) + "\t"
                                  + std::to_string(fp) + "\t" + std::to_string(fn) + "\t");
            outStep += (std::to_string(iters) + "\t" + std::to_string(time(nullptr) - start) + "\t" +
                        std::to_string(cutoff) + "\t" + std::to_string(numBins) + "\t" + std::to_string(cutoff) + "\t")
                    + std::to_string(tp) + '\t' + std::to_string(tn) + '\t' + std::to_string(fp) + '\t' +
                    std::to_string(fn) +
                    '\t';
            for (double result : results) {
                clusterOutputData += (std::to_string(result) + "\t");
                outStep += std::to_string(result) + "\t";
            }
            outStep += "\n";
        }
        ListVector *list = nullptr;

        clusterOutputData += "\n\n";
        list = cluster.getList();
        list->setLabel(std::to_string(cutoff));
        //
        if (printHeaders) {
            //only print headers the first time
            printHeaders = false;
        } else { list->setPrintedLabels(printHeaders); }
        clusterMatrixOutput = list->print(listFile);
#if !DEBUG_RCPP
        std::ofstream outPutWriter("F:\\muMS2_Clustering\\muMS2_Clustering\\Distance_files\\cluster_ouput.list");
        while(!outPutWriter.is_open()) {
            TestHelper::Print("Enter new path: ");
            std::string filePath;
            std::cin >> filePath;
            outPutWriter.open(filePath);
        }
        outPutWriter << clusterMatrixOutput;
        outPutWriter.close();
#endif
        delete list;
        results = cluster.getStats(tp, tn, fp, fn);

        sensFile += std::to_string(cutoff) + '\t' + std::to_string(cutoff) + '\t' + std::to_string(tp) + '\t' +
                std::to_string(tn) + '\t' +
                std::to_string(fp) + '\t' + std::to_string(fn) + '\t';
        for (double result : results) { sensFile + std::to_string(result) + '\t'; }
    }
    delete matrix;
    return 0;
}
