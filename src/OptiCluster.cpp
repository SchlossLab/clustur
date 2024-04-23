//
// Created by Gregory Johnson on 3/29/24.
//

#include "MothurDependencies/OptiCluster.h"
#include "MothurDependencies/Metrics/accuracy.hpp"
#include "MothurDependencies/Metrics/f1score.hpp"
#include "MothurDependencies/Metrics/fdr.hpp"
#include "MothurDependencies/Metrics/mcc.hpp"
#include "MothurDependencies/Metrics/npv.hpp"
#include "MothurDependencies/Metrics/ppv.hpp"
#include "MothurDependencies/Metrics/sensitivity.hpp"
#include "MothurDependencies/Metrics/specificity.hpp"

OptiCluster::OptiCluster(OptiData *mt, ClusterMetric *met, long long ns) {
    matrix = mt;
    metric = met;
    truePositives = 0;
    trueNegatives = 0;
    falseNegatives = 0;
    falsePositives = 0;
    numSingletons = ns;
}

/***********************************************************************/
//randomly assign sequences to OTUs
int OptiCluster::initialize(double &value, bool randomize, std::string initialize) {
    numSeqs = matrix->getNumSeqs();
    truePositives = 0;
    falsePositives = 0;
    falseNegatives = 0;
    trueNegatives = 0;

    bins.resize(numSeqs); //place seqs in own bin

    std::vector<long long> temp;
    bins.push_back(temp);
    seqBin[numSeqs] = -1;
    insertLocation = numSeqs;
    Utils util;

    if (initialize == "singleton") {
        //put everyone in own bin
        for (int i = 0; i < numSeqs; i++) { bins[i].push_back(i); }

        //maps randomized sequences to bins
        for (int i = 0; i < numSeqs; i++) {
            seqBin[i] = bins[i][0];
            randomizeSeqs.push_back(i);
        }

        if (randomize) { util.mothurRandomShuffle(randomizeSeqs); }

        //for each sequence (singletons removed on read)
        for (std::map<long long, long long>::iterator it = seqBin.begin(); it != seqBin.end(); it++) {
            if (it->second == -1) {
            } else {
                long long numCloseSeqs = (matrix->getNumClose(it->first)); //does not include self
                falseNegatives += numCloseSeqs;
            }
        }
        falseNegatives /= 2; //square matrix
        trueNegatives = numSeqs * (numSeqs - 1) / 2 - (falsePositives + falseNegatives + truePositives);
        //since everyone is a singleton no one clusters together. True negative = num far apart
    } else {
        //put everyone in first bin
        for (int i = 0; i < numSeqs; i++) {
            bins[0].push_back(i);
            seqBin[i] = 0;
            randomizeSeqs.push_back(i);
        }

        if (randomize) { util.mothurRandomShuffle(randomizeSeqs); }

        //for each sequence (singletons removed on read)
        for (std::map<long long, long long>::iterator it = seqBin.begin(); it != seqBin.end(); it++) {
            if (it->second == -1) {
            } else {
                long long numCloseSeqs = (matrix->getNumClose(it->first)); //does not include self
                truePositives += numCloseSeqs;
            }
        }
        truePositives /= 2; //square matrix
        falsePositives = numSeqs * (numSeqs - 1) / 2 - (trueNegatives + falseNegatives + truePositives);
    }

    value = metric->getValue(truePositives, trueNegatives, falsePositives, falseNegatives);

    return value;
}

/***********************************************************************/
/* for each sequence with mutual information (close)
 * remove from current OTU and calculate MCC when sequence forms its own OTU or joins one of the other OTUs where there is a sequence within the `threshold` (no need to calculate MCC if the paired sequence is already in same OTU and no need to try every OTU - just those where there's a close sequence)
 * keep or move the sequence to the OTU where the `metric` is the largest - flip a coin on ties */
bool OptiCluster::update(double &listMetric) {
    //for each sequence (singletons removed on read)
    for (int i = 0; i < randomizeSeqs.size(); i++) {
        const std::map<long long, long long>::iterator it = seqBin.find(randomizeSeqs[i]);

        long long seqNumber = it->first;
        long long binNumber = it->second;

        if (binNumber == -1) {
        } else {
            double tn, tp, fp, fn;
            double bestMetric = -1;
            double bestBin, bestTp, bestTn, bestFn, bestFp;
            tn = trueNegatives;
            tp = truePositives;
            fp = falsePositives;
            fn = falseNegatives;

            //close / far count in current bin
            std::vector<double> results = getCloseFarCounts(seqNumber, binNumber);
            double cCount = results[0];
            double fCount = results[1];

            //metric in current bin
            bestMetric = metric->getValue(tp, tn, fp, fn);
            bestBin = binNumber;
            bestTp = tp;
            bestTn = tn;
            bestFp = fp;
            bestFn = fn;

            //if not already singleton, then calc value if singleton was created
            if (!((bins[binNumber].size()) == 1)) {
                //make a singleton
                //move out of old bin
                fn += cCount;
                tn += fCount;
                fp -= fCount;
                tp -= cCount;
                double singleMetric = metric->getValue(tp, tn, fp, fn);
                if (singleMetric > bestMetric) {
                    bestBin = -1;
                    bestTp = tp;
                    bestTn = tn;
                    bestFp = fp;
                    bestFn = fn;
                    bestMetric = singleMetric;
                }
            }

            std::set<long long> binsToTry;
            std::set<long long> closeSeqs = matrix->getCloseSeqs(seqNumber);
            for (std::set<long long>::iterator itClose = closeSeqs.begin(); itClose != closeSeqs.end(); itClose++) {
                binsToTry.insert(seqBin[*itClose]);
            }

            //merge into each "close" otu
            for (std::set<long long>::iterator it = binsToTry.begin(); it != binsToTry.end(); it++) {
                tn = trueNegatives;
                tp = truePositives;
                fp = falsePositives;
                fn = falseNegatives;
                fn += cCount;
                tn += fCount;
                fp -= fCount;
                tp -= cCount; //move out of old bin
                results = getCloseFarCounts(seqNumber, *it);
                fn -= results[0];
                tn -= results[1];
                tp += results[0];
                fp += results[1]; //move into new bin
                double newMetric = metric->getValue(tp, tn, fp, fn); //score when sequence is moved
                //new best
                if (newMetric > bestMetric) {
                    bestMetric = newMetric;
                    bestBin = (*it);
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
    }

    listMetric = metric->getValue(truePositives, trueNegatives, falsePositives, falseNegatives);


    return 0;
}

/***********************************************************************/
std::vector<double> OptiCluster::getCloseFarCounts(long long seq, long long newBin) {
    std::vector<double> results;
    results.push_back(0);
    results.push_back(0);
    if (newBin == -1) {
    } //making a singleton bin. Close but we are forcing apart.
    else {
        //merging a bin
        for (int i = 0; i < bins[newBin].size(); i++) {
            if (seq == bins[newBin][i]) {
            } //ignore self
            else if (!matrix->isClose(seq, bins[newBin][i])) { results[1]++; }
            //this sequence is "far away" from sequence i - above the cutoff
            else { results[0]++; } //this sequence is "close" to sequence i - distance between them is less than cutoff
        }
    }
    return results;
}

//TODO figure out why the bad-allocation execption is being thrown
ListVector *OptiCluster::getList() {
    ListVector *list = new ListVector();
    ListVector *singleton = matrix->getListSingle();
    // TestHelper::Print("Made it here, listVector\n");
    if (singleton != nullptr) {
        //add in any sequences above cutoff in read. Removing these saves clustering time.
        for (int i = 0; i < singleton->getNumBins(); i++) {
            if (singleton->get(i) != "") {
                list->push_back(singleton->get(i));
            }
        }
        delete singleton;
    }
    // TestHelper::Print("Made it here, listVector2\n");
    const std::string sizeString = "Size: " + std::to_string(bins.size());
    for (int i = 0; i < bins.size(); i++) {
        if (bins[i].size() != 0) {
            std::string otu = matrix->getName(bins[i][0]);
            //   TestHelper::Print(otu + "\n");
            for (int j = 1; j < bins[i].size(); j++) {
                // TestHelper::Print(i + "\n");
                otu += "," + matrix->getName(bins[i][j]);
            }
            list->push_back(otu);
        }
    }
    return list;
}

/***********************************************************************/
std::vector<double> OptiCluster::getStats(double &tp, double &tn, double &fp, double &fn) {
    double singletn = matrix->getNumSingletons() + numSingletons;
    double tempnumSeqs = numSeqs + singletn;

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

std::vector<double> OptiCluster::getCloseFarFitCounts(long long seq, long long newBin) {
    std::vector<double> results;
    results.push_back(0);
    results.push_back(0); //results[0] = close count, results[1] = far count

    if (newBin == -1) {
    } //making a singleton bin. Close but we are forcing apart.
    else {
        //merging a bin
        for (long long i = 0; i < bins[newBin].size(); i++) {
            if (seq == bins[newBin][i]) {
            } //ignore self
            else {
                bool isFit = true;
                bool closeFit = matrix->isCloseFit(seq, bins[newBin][i], isFit);
                if (closeFit) {
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
long long OptiCluster::getNumBins() {
    long long singletn = matrix->getNumSingletons();

    for (int i = 0; i < bins.size(); i++) {
        if (bins[i].size() != 0) {
            singletn++;
        }
    }

    return singletn;
}

/***********************************************************************/
long long OptiCluster::findInsert() {
    //initially there are bins for each sequence (excluding singletons removed on read)
    for (long long i = 0; i < bins.size(); i++) {
        if (bins[i].size() == 0) { return i; } //this bin is empty
    }

    return -1;
}

/***********************************************************************/
