//
// Created by Gregory Johnson on 3/29/24.
//
#include <Rcpp.h>
#include <unordered_set>
#include <fstream>
#include "Clusters/OptiCluster.h"
#include "Clusters/Metrics/accuracy.h"
#include "Clusters/Metrics/f1score.h"
#include "Clusters/Metrics/fdr.h"
#include "Clusters/Metrics/mcc.h"
#include "Clusters/Metrics/npv.h"
#include "Clusters/Metrics/ppv.h"
#include "Clusters/Metrics/sensitivity.h"
#include "Clusters/Metrics/specificity.h"
#include "DataExporters/OptiClusterData.h"
#include "Clusters/Metrics/tptn.h"
#include <chrono>



OptiCluster::OptiCluster(OptiData *mt, ClusterMetric *met, const double cutoff, const long long ns) : matrix(mt),
    metric(met), numSeqs(0), insertLocation(0), numSingletons(ns), fittruePositives(0), fittrueNegatives(0),
    fitfalsePositives(0),
    fitfalseNegatives(0),
    combotruePositives(0),
    combotrueNegatives(0),
    combofalsePositives(0),
    combofalseNegatives(0),
    numFitSeqs(0), numFitSingletons(0),
    numComboSeqs(0),
    numComboSingletons(0),
    cutoff(cutoff) {
    truePositives = 0;
    trueNegatives = 0;
    falseNegatives = 0;
    falsePositives = 0;
}

OptiCluster::~OptiCluster() {
    delete metric;
    delete matrix;
}

/***********************************************************************/
//randomly assign sequences to OTUs
int OptiCluster::initialize(double &value, const bool randomize, const std::string& initialize) {
    numSeqs = matrix->getNumSeqs();
    truePositives = 0;
    falsePositives = 0;
    falseNegatives = 0;
    trueNegatives = 0;

    bins.resize(numSeqs); //place seqs in own bin

    const std::vector<long long> temp;
    seqBin = std::vector<long long>(numSeqs, 0);
    bins.push_back(temp);
    // seqBin[numSeqs] = -1;
    insertLocation = numSeqs;

    if (initialize == "singleton") {
        //put everyone in own bin
        for (int i = 0; i < numSeqs; i++) {
            bins[i].push_back(i);
            seqBin[i] = i;
            randomizeSeqs.push_back(i);
        }

        //maps randomized sequences to bins
        // for (int i = 0; i < numSeqs; i++) {
        //     seqBin[i] = bins[i][0];
        //     // index[
        //     randomizeSeqs.push_back(i);
        // }

        if (randomize) { Utils::mothurRandomShuffle(randomizeSeqs); }

        //for each sequence (singletons removed on read)
        for (const auto seq: seqBin) {
            const long long numCloseSeqs = (matrix->getNumClose(seq)); //does not include self
            falseNegatives += static_cast<double>(numCloseSeqs);
        }
        falseNegatives /= 2; //square matrix
        const auto nSeqs = static_cast<double>(numSeqs);
        trueNegatives = nSeqs * (nSeqs - 1) / 2 - (falsePositives + falseNegatives + truePositives);
        //since everyone is a singleton no one clusters together. True negative = num far apart
    } else {
        //put everyone in first bin
        for (int i = 0; i < numSeqs; i++) {
            bins[0].push_back(i);
            seqBin[i] = 0;
            randomizeSeqs.push_back(i);
        }

        if (randomize) { Utils::mothurRandomShuffle(randomizeSeqs); }

        //for each sequence (singletons removed on read)
        for (const auto seq : seqBin) {
            const long long numCloseSeqs = (matrix->getNumClose(seq)); //does not include self
            truePositives += static_cast<double>(numCloseSeqs);
        }
        truePositives /= 2; //square matrix
        const auto nSeqs = static_cast<double>(numSeqs);
        falsePositives = nSeqs * (nSeqs - 1) / 2 - (trueNegatives + falseNegatives + truePositives);
    }

    value = metric->getValue(truePositives, trueNegatives, falsePositives, falseNegatives);

    return static_cast<int>(value);
}

/***********************************************************************/
/* for each sequence with mutual information (close)
 * remove from current OTU and calculate MCC when sequence forms its own OTU or joins one of the other OTUs where there is a sequence within the `threshold` (no need to calculate MCC if the paired sequence is already in same OTU and no need to try every OTU - just those where there's a close sequence)
 * keep or move the sequence to the OTU where the `metric` is the largest - flip a coin on ties */
bool OptiCluster::update(double &listMetric) {
    //for each sequence (singletons removed on read)
    for (const int randomizeSeq : randomizeSeqs) {
        // const auto it = seqBin.find(randomizeSeq);
        // seqBin[randomizeSeq]

        long long seqNumber = randomizeSeq;
        const std::vector<long long> closeSeqs = matrix->getCloseSeqs(seqNumber);
        const long long binNumber = seqBin[randomizeSeq];

        double tn = trueNegatives;
        double tp = truePositives;
        double fp = falsePositives;
        double fn = falseNegatives;

        //close / far count in current bin
        std::vector<double> results = getCloseFarCounts(seqNumber, binNumber);
        const double cCount = results[0];
        const double fCount = results[1];

        //metric in current bin
        double bestMetric = metric->getValue(tp, tn, fp, fn);
        long long bestBin = binNumber;
        double bestTp = tp;
        double bestTn = tn;
        double bestFp = fp;
        double bestFn = fn;

        //if not already singleton, then calc value if singleton was created
        if ((bins[binNumber].size()) != 1) {
            //make a singleton
            //move out of old bin
            fn += cCount;
            tn += fCount;
            fp -= fCount;
            tp -= cCount;
            if (const double singleMetric = metric->getValue(tp, tn, fp, fn); singleMetric > bestMetric) {
                bestBin = -1;
                bestTp = tp;
                bestTn = tn;
                bestFp = fp;
                bestFn = fn;
                bestMetric = singleMetric;
            }
        }



        std::vector<long long> binsToTry;
        //binsToTry.reserve(closeSeqs.size());
        for (const long long closeSeq : closeSeqs) {
            binsToTry.emplace_back(seqBin[closeSeq]);
        }
        std::sort(binsToTry.begin(), binsToTry.end());
        const auto it = std::unique(binsToTry.begin(), binsToTry.end());
        binsToTry.resize(std::distance(binsToTry.begin(), it));

        //merge into each "close" otu
        for (const long long bin : binsToTry) {
            tn = trueNegatives;
            tp = truePositives;
            fp = falsePositives;
            fn = falseNegatives;
            fn += cCount;
            tn += fCount;
            fp -= fCount;
            tp -= cCount; //move out of old bin
            results = getCloseFarCounts(seqNumber, bin);
            fn -= results[0];
            tn -= results[1];
            tp += results[0];
            fp += results[1]; //move into new bin
            //new best
            if (const double newMetric = metric->getValue(tp, tn, fp, fn); newMetric > bestMetric) {
                bestMetric = newMetric;
                bestBin = bin;
                bestTp = tp;
                bestTn = tn;
                bestFp = fp;
                bestFn = fn;
            }
        }

        bool usedInsert = false;
        if (bestBin == -1) {
            bestBin = insertLocation;
            usedInsert = true;
        }

        if (bestBin != binNumber) {
            truePositives = bestTp;
            trueNegatives = bestTn;
            falsePositives = bestFp;
            falseNegatives = bestFn;

            //move seq from i to j
            bins[bestBin].push_back(seqNumber); //add seq to bestbin
            bins[binNumber].erase(remove(bins[binNumber].begin(), bins[binNumber].end(), seqNumber),
                                  bins[binNumber].end()); //remove from old bin i
        }

        if (usedInsert) { insertLocation = findInsert(); }

        //update seqBins
        seqBin[seqNumber] = bestBin; //set new OTU location
    }


    listMetric = metric->getValue(truePositives, trueNegatives, falsePositives, falseNegatives);


    return false;
}

/***********************************************************************/
std::vector<double> OptiCluster::getCloseFarCounts(const long long seq, const long long newBin) const {
    std::vector<double> results(2, 0);
    if (newBin == -1) {
    } //making a singleton bin. Close but we are forcing apart.
    else {
        //merging a bin
        for (const long long bin : bins[newBin]) {
            if (seq == bin) continue; //ignore self
            results[!matrix->isClose(seq, bin)]++;
            // else if (!matrix->isClose(seq, bin)) { results[1]++; }
            // //this sequence is "far away" from sequence i - above the cutoff
            // else { results[0]++; } //this sequence is "close" to sequence i - distance between them is less than cutoff
        }
    }
    return results;
}

//TODO figure out why the bad-allocation execption is being thrown
ListVector OptiCluster::getList() const {
    ListVector list;
    const ListVector *singleton = matrix->getListSingle();
    // TestHelper::Print("Made it here, listVector\n");
    if (singleton != nullptr) {
        //add in any sequences above cutoff in read. Removing these saves clustering time.
        for (int i = 0; i < singleton->getNumBins(); i++) {
            if (!singleton->get(i).empty()) {
                list.push_back(singleton->get(i));
            }
        }
        delete singleton;
    }
    // TestHelper::Print("Made it here, listVector2\n");
    const std::string sizeString = "Size: " + std::to_string(bins.size());
    for (const auto& bin : bins){
        if (!bin.empty()) {
            std::string otu = matrix->getName(bin[0]);
            //   TestHelper::Print(otu + "\n");
            for (size_t j = 1; j < bin.size(); j++) {
                // TestHelper::Print(i + "\n");
                otu += "," + matrix->getName(bin[j]);
            }
            list.push_back(otu);
        }
    }
    return list;
}

/***********************************************************************/
std::vector<double> OptiCluster::getStats(double &tp, double &tn, double &fp, double &fn) const {
    long long singletn = matrix->getNumSingletons() + numSingletons;
    auto tempnumSeqs = static_cast<double>(numSeqs + singletn);

    tp = truePositives;
    fp = falsePositives;
    fn = falseNegatives;
    tn = tempnumSeqs * (tempnumSeqs - 1) / 2 - (falsePositives + falseNegatives + truePositives);
    //adds singletons to tn

    std::vector<double> results;

    Sensitivity sens;
    double sensitivity = sens.getValue(tp, tn, fp, fn);
    results.push_back(sensitivity);
    Specificity spec;
    double specificity = spec.getValue(tp, tn, fp, fn);
    results.push_back(specificity);
    PPV ppv;
    double positivePredictiveValue = ppv.getValue(tp, tn, fp, fn);
    results.push_back(positivePredictiveValue);
    NPV npv;
    double negativePredictiveValue = npv.getValue(tp, tn, fp, fn);
    results.push_back(negativePredictiveValue);
    FDR fdr;
    double falseDiscoveryRate = fdr.getValue(tp, tn, fp, fn);
    results.push_back(falseDiscoveryRate);
    Accuracy acc;
    double accuracy = acc.getValue(tp, tn, fp, fn);
    results.push_back(accuracy);
    MCC mcc;
    double matthewsCorrCoef = mcc.getValue(tp, tn, fp, fn);
    results.push_back(matthewsCorrCoef);
    F1Score f1;
    double f1Score = f1.getValue(tp, tn, fp, fn);
    results.push_back(f1Score);

    return results;
}

std::vector<double> OptiCluster::getCloseFarFitCounts(const long long seq, const long long newBin) const {
    std::vector<double> results;
    results.push_back(0);
    results.push_back(0); //results[0] = close count, results[1] = far count

    if (newBin == -1) {
    } //making a singleton bin. Close but we are forcing apart.
    else {
        //merging a bin
        for (const long long& bin : bins[newBin]) {
            if (seq == bin) {
            } //ignore self
            else {
                bool isFit = true;
                if (const bool closeFit = matrix->isCloseFit(seq, bin, isFit)) {
                    //you are close if you are fit and close
                    results[0]++;
                } else if (isFit) { results[1]++; } //this sequence is "far away" and fit - above the cutoff
            }
        }
    }

    return results;
}

/***********************************************************************/

/***********************************************************************/
long long OptiCluster::getNumBins() const {
    long long singleton = matrix->getNumSingletons();

    for (const auto & bin : bins) {
        if (!bin.empty()) {
            singleton++;
        }
    }

    return singleton;
}

/***********************************************************************/
long long OptiCluster::findInsert() const {
    //initially there are bins for each sequence (excluding singletons removed on read)
    for (size_t i = 0; i < bins.size(); i++) {
        if (bins[i].empty()) { return static_cast<long long>(i); } //this bin is empty
    }

    return -1;
}

/***********************************************************************/



ClusterExport* OptiCluster::Execute() {
    // clusterMetrics += ("\nClustering " + distfile + "\n");
    std::string initializeString = "singleton";
    auto* data = new OptiClusterData("");
    std::string clusterMetrics;
    std::string sensFile;
    bool canShuffle = true;
    double stableMetric = 0;
    int maxIters = 100;
    std::ofstream listFile;
    std::vector<std::string> sensfileHeaders{"label","cutoff","ttp","tn","fp","fn","sensitivity",
        "specificity","ppv","npv","fdr","accuracy","mcc","f1score"};
    // sensFile = "label\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n";

    std::vector<std::string> clusterMetricsHeaders{"iter","time","label","num_otus","cutoff","tp","tn",
      "fp","fn","sensitivity","specificity","ppv","npv",
        "fdr", "accuracy", "mcc", "f1score"};
    // clusterMetrics = (
    //     "iter\ttime\tlabel\tnum_otus\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n");
    bool printHeaders = true;

    if (!matrix->mccValidCalc()) {
        Rcpp::warning("[WARNING]: The mcc metric is not suitible for your data with a cutoff of " +
            std::to_string(cutoff) + " using tptn instead.");
        delete metric;
        metric = new TPTN();
    }

    int iters = 0;
    double listVectorMetric = 0; //worst state
    double delta = 1;
    long long numBins;
    double tp, tn, fp, fn;
    std::vector<double> stats;
    std::vector<std::string> clusterMetricList;
    initialize(listVectorMetric, canShuffle, initializeString);
    stats = getStats(tp, tn, fp, fn);
    numBins = getNumBins();
    clusterMetrics = ("0,0," + std::to_string(cutoff) + "," + std::to_string(numBins) + "," +
                       std::to_string(cutoff) + "," + std::to_string(tp) + "," + std::to_string(tn) + "," +
                       std::to_string(fp) + "," + std::to_string(fn) + ",");


    for (double result: stats) {
        clusterMetrics += (std::to_string(result) + ",");
    }
    Utils::AddRowToDataFrameMap(dataframeMapClusterMetrics, clusterMetrics, clusterMetricsHeaders);

    while ((delta > stableMetric) && (iters < maxIters)) {
        //long start = std::time(nullptr);
        double oldMetric = listVectorMetric;
        auto startTime = std::chrono::system_clock::now();
        update(listVectorMetric);

        delta = std::abs(oldMetric - listVectorMetric);
        iters++;

        stats = getStats(tp, tn, fp, fn);

        numBins = getNumBins();
        auto endTime = std::chrono::system_clock::now();
        std::chrono::duration<double> currentTime = endTime - startTime;
        clusterMetrics = (std::to_string(iters) + "," + std::to_string(currentTime.count()) + "," +
                           std::to_string(cutoff) + "," + std::to_string(numBins) + "," +
                           std::to_string(cutoff) + "," + std::to_string(tp) + "," + std::to_string(tn) + ","
                           + std::to_string(fp) + "," + std::to_string(fn) + ",");
        for (double result: stats) {
            clusterMetrics += (std::to_string(result) + ",");
        }
        Utils::AddRowToDataFrameMap(dataframeMapClusterMetrics, clusterMetrics, clusterMetricsHeaders);
    }

    ListVector list = getList();
    list.setLabel(std::to_string(cutoff));
    OptiClusterInformation clusterInformation;
    clusterInformation.label = std::to_string(cutoff);
    clusterInformation.numberOfOtu = static_cast<int>(numBins);
    clusterInformation.clusterBins = list.print(listFile);
    data->AddToData(clusterInformation);
    data->SetListVector(list, std::to_string(cutoff));
    stats = getStats(tp, tn, fp, fn);

    sensFile += std::to_string(cutoff) + ',' + std::to_string(cutoff) + ',' + std::to_string(tp) + ',' +
            std::to_string(tn) + ',' +
            std::to_string(fp) + ',' + std::to_string(fn) + ',';
    for (double result: stats) { sensFile += std::to_string(result) + ','; }
    Utils::AddRowToDataFrameMap(dataframeMapSensMetrics, sensFile, sensfileHeaders);

    return data;
}