

//
//  optifitcluster.cpp
//  Mothur
//
//  Created by Sarah Westcott on 5/10/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//


#include "Clusters/Optifitcluster.h"
#include "Clusters/OptiCluster.h"
#include "Clusters/Metrics/accuracy.h"
#include "Clusters/Metrics/f1score.h"
#include "Clusters/Metrics/fdr.h"
#include "Clusters/Metrics/mcc.h"
#include "Clusters/Metrics/npv.h"
#include "Clusters/Metrics/ppv.h"
#include "Clusters/Metrics/sensitivity.h"
#include "Clusters/Metrics/specificity.h"
#include "DataExporters/OptifitClusterData.h"
#include "DataStructures/OptiRefMatrix.h"
#include <chrono>

#include "Clusters/Cluster.h"


/***********************************************************************/
OptiFitCluster::OptiFitCluster(OptiRefMatrix* mt, ClusterMetric* met, const double cutoff, const long long ns)
    : metric(met), matrix(mt), cutoff(cutoff), numComboSingletons(ns) {
    maxRefBinNumber = 0;
    closed = false;

    numFitSeqs = 0;  fittruePositives = 0; fitfalsePositives = 0; fitfalseNegatives = 0; fittrueNegatives = 0; numFitSingletons = 0;
    numComboSeqs = 0; numComboSingletons = 0; combotruePositives = 0; combofalsePositives = 0; combofalseNegatives = 0; combotrueNegatives = 0;
}

OptiFitCluster::OptiFitCluster(OptiRefMatrix* mt, ClusterMetric* met, const ListVector& listVector, const double cutoff, const long long ns)
    : metric(met), matrix(mt), listVector(listVector), cutoff(cutoff), numComboSingletons(ns) {
    maxRefBinNumber = 0;
    closed = false;

    numFitSeqs = 0;  fittruePositives = 0; fitfalsePositives = 0; fitfalseNegatives = 0; fittrueNegatives = 0; numFitSingletons = 0;
    numComboSeqs = 0; numComboSingletons = 0; combotruePositives = 0; combofalsePositives = 0; combofalseNegatives = 0; combotrueNegatives = 0;
}

void OptiFitCluster::Reset() {
    maxRefBinNumber = 0;
    closed = false;

    numFitSeqs = 0;  fittruePositives = 0; fitfalsePositives = 0; fitfalseNegatives = 0; fittrueNegatives = 0; numFitSingletons = 0;
    numComboSeqs = 0; numComboSingletons = 0; combotruePositives = 0; combofalsePositives = 0; combofalseNegatives = 0; combotrueNegatives = 0;
}


ClusterExport * OptiFitCluster::Execute() {
    time_t estart = time(nullptr);
    constexpr bool selfReference = false;
    // if (metricName == "mcc")             { metric = new MCC();              }
    // else if (metricName == "sens")       { metric = new Sensitivity();      }
    // else if (metricName == "spec")       { metric = new Specificity();      }
    // else if (metricName == "tptn")       { metric = new TPTN();             }
    // else if (metricName == "tp")         { metric = new TP();               }
    // else if (metricName == "tn")         { metric = new TN();               }
    // else if (metricName == "fp")         { metric = new FP();               }
    // else if (metricName == "fn")         { metric = new FN();               }
    // else if (metricName == "f1score")    { metric = new F1Score();          }
    // else if (metricName == "accuracy")   { metric = new Accuracy();         }
    // else if (metricName == "ppv")        { metric = new PPV();              }
    // else if (metricName == "npv")        { metric = new NPV();              }
    // else if (metricName == "fdr")        { metric = new FDR();              }
    // else if (metricName == "fpfn")       { metric = new FPFN();             }

    std::map<std::string, int> counts;
    // std::string dupsFile = countfile; nameOrCount = "count";
    // if (namefile != "") { dupsFile = namefile; nameOrCount = "name"; }
    // else { CountTable ct; ct.readTable(countfile, false, false); counts = ct.getNameMap();  }

    // if (outputdir == "") { outputdir += util.hasPath(distfile); }
    // fileroot = outputdir + util.getRootName(util.getSimpleName(distfile));

    std::string listFile = ""; std::string bestListFileName = ""; std::string outputName = "";

    if (selfReference) { //de novo
        // std::map<std::string, std::string> variables;
        // variables["[filename]"] = fileroot;
        // variables["[clustertag]"] = "optifit_" + metric->getName();
        // outputName = getOutputFileName("steps", variables);
        //
        // if ((accnosfile == "") && (!createAccnos)) { //denovo with mothur randomly assigning references
        //
        // Rcpp::message("\nRandomly assigning reads from the reference  as reference sequences\n");

        //distfile, distFormat, dupsFile, dupsFormat, cutoff, percentage to be fitseqs - will randomly assign as fit
        // OptiData* matrix = new OptiRefMatrix(distfile, "column", dupsFile, nameOrCount, cutoff, fitPercent, refWeight);

        return runDenovoOptiCluster(counts, outputName);

        // std::string sensspecFilename = fileroot+ tag + ".sensspec";
        // ofstream sensFile;
        // util.openOutputFile(sensspecFilename,    sensFile);
        // outputNames.push_back(sensspecFilename); outputTypes["sensspec"].push_back(sensspecFilename);

        //evaluate results
        // bestListFileName = compareSensSpec(matrix, metric, sensFile);
        //
        //     delete matrix;

    }else { //reference with accnos file or reference list file assigning references

        std::unordered_set<std::string> refNames; std::vector<std::string> refLabels; std::vector< std::vector<std::string> > otus;

        // if (accnosfile != "") { //use accnos file to assign references
        //
        //     m->mothurOut("\nUsing sequences from " + accnosfile + " as reference sequences\n");
        //
        //     refNames = Utils::readAccnos(accnosfile);

        // }else if (createAccnos) { //assign references based on reflist parameter

        // m->mothurOut("\nUsing OTUs from " + reflistfile + " as reference OTUs\n");
        //
        // InputData input(reflistfile, "list", nullVector);
        // std::set<std::string> processedLabels, userLabels;
        // std::string lastLabel = "";
        //
        // ListVector* reflist = Utils::getNextList(input, true, userLabels, processedLabels, lastLabel);
        //
        // refLabels = reflist->getLabels();
        // for (int i = 0; i < refLabels.size(); i++) { refLabels[i] = "Ref_" + refLabels[i];  }
        //
        // refNames = util.getSetFromList(reflist, otus); delete reflist;
        //   }

        //distfile, distFormat, dupsFile, dupsFormat, cutoff, accnos containing refseq name
        // OptiData* matrix = new OptiRefMatrix(distfile, "column", dupsFile, nameOrCount, cutoff, refNames);

        //fit seqs
        // return runUserRefOptiCluster(metric,refLabels, otus);
        //
        // ofstream listFile; string listFileName = fileroot+ tag + ".list";
        // util.openOutputFile(listFileName,    listFile);
        //
        // if(countfile != "") { list->print(listFile, counts); }
        // else { list->print(listFile); }
        // listFile.close();
        //
        // listFiles.push_back(listFileName);
        // bestListFileName = listFileName;
        //
        // delete list;

    }//reference with files containing reference seqs

    // createReferenceNameCount(); //creates reference name or count file if needed
    //
    // std::string distanceFile = calcDists();  //calc distance matrix for fasta file and distances between fasta file and reffasta file
    //
    // // if (outputdir == "") { outputdir += util.hasPath(distanceFile); }
    // // fileroot = outputdir + util.getRootName(util.getSimpleName(distanceFile));
    //
    // std::map<std::string, std::string> variables;
    // variables["[filename]"] = fileroot;
    // variables["[clustertag]"] = "optifit_" + metric->getName();
    // // outputName = getOutputFileName("steps", variables);
    //
    // // m->mothurOut("\nUsing OTUs from " + reflistfile + " as reference OTUs\n");
    //
    // //calc sens.spec values for reference
    // // InputData input(reflistfile, "list", nullVector);
    // ListVector* list = input.getListVector();
    //
    // //add tag to OTULabels to indicate the reference
    // std::vector<std::string> refListLabels = list->getLabels();
    // for (int i = 0; i < refListLabels.size(); i++) { refListLabels[i] = "Ref_" + refListLabels[i];  }
    // list->setLabels(refListLabels);
    //
    // string refDupsFile = refcountfile;
    // if (refNameOrCount == "name") { refDupsFile = refnamefile; }
    //
    // OptiData* matrix = new OptiRefMatrix(refdistfile, refDupsFile, refNameOrCount, refformat, cutoff, distfile, dupsFile, nameOrCount, "column", comboDistFile, "column");
    //

    return runRefOptiCluster(listVector, counts, outputName);
    // // listFiles.push_back(listFile);
    //
    // // bestListFileName = listFile;
    // //
    // delete matrix;

    // delete metric;

    // if (m->getControl_pressed()) {     for (int j = 0; j < outputNames.size(); j++) { util.mothurRemove(outputNames[j]); }  return 0; }
    //
    // outputNames.push_back(outputName); outputTypes["steps"].push_back(outputName);
    // outputNames.push_back(bestListFileName); outputTypes["list"].push_back(bestListFileName);
    //
    // if (m->getControl_pressed()) {     for (int j = 0; j < outputNames.size(); j++) { util.mothurRemove(outputNames[j]); }  return 0; }
    //
    // m->mothurOut("It took " + toString(time(nullptr) - estart) + " seconds to fit sequences to reference OTUs.\n");
    //
    // //set list file as new current listfile
    // string currentName = "";
    // itTypes = outputTypes.find("list");
    // if (itTypes != outputTypes.end()) {
    //     if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setListFile(currentName); }
    // }
    //
    // itTypes = outputTypes.find("accnos");
    // if (itTypes != outputTypes.end()) {
    //     if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setAccnosFile(currentName); }
    // }
    //
    // m->mothurOut("\nOutput File Names: \n");
    // for (int i = 0; i < outputNames.size(); i++) {    m->mothurOut(outputNames[i]+"\n");     }
    // m->mothurOutEndLine();

    return nullptr;
}

/***********************************************************************/
int OptiFitCluster::initialize(double& value, const bool randomize, std::vector<std::vector<std::string > >& existingBins,
                               const std::vector<std::string>& bls, const std::string& meth, const bool denov) {
    double reftruePositives, reftrueNegatives, reffalsePositives, reffalseNegatives, numRefSeqs;
    numRefSeqs = 0; reftruePositives = 0; reffalsePositives = 0; reffalseNegatives = 0; reftrueNegatives = 0;

    if (meth == "closed") { closed = true; }
    denovo = denov;

    std::vector< std::vector< long long> > translatedBins;
    randomizeSeqs = matrix->getTranslatedBins(existingBins, translatedBins); //otus in existingBins, otus with matrix names
    // Bins and translated bins are the same thing...
    int binNumber = 0;
    int placeHolderIndex = -1;
    bins.clear();
    bins.reserve(translatedBins.size() + randomizeSeqs.size());
    seqBin.clear();
    for (long long i = 0; i < translatedBins.size(); i++) {
        binLabels[binNumber] = bls[i];
        // bins[i] = translatedBins[i];
        bins.emplace_back(translatedBins[i]);
        numRefSeqs += translatedBins[i].size();

        for (int j = 0; j < translatedBins[i].size(); j++) {
            for (int k = 0; k < j; k++) {
                if (translatedBins[i][j] < 0) { //no dists in matrix
                    translatedBins[i][j] = placeHolderIndex; placeHolderIndex--;
                    reffalsePositives++;
                }else { //j has distances in the matrix, but is it close to k?
                    if (matrix->isClose(translatedBins[i][j], translatedBins[i][k])) {
                        reftruePositives++;
                    }else { reffalsePositives++; }
                }
            }
            seqBin[translatedBins[i][j]] = binNumber;
        }
        binNumber++;
    }

    maxRefBinNumber = binNumber;
    reffalseNegatives = matrix->getNumRefDists() - reftruePositives; //number of distance in matrix for reference seqs - reftruePositives
    reftrueNegatives = numRefSeqs * (numRefSeqs-1)/2 - (reffalsePositives + reffalseNegatives + reftruePositives);

    //add fit seqs as singletons
    const int numRefBins = translatedBins.size();
    numFitSingletons = 0;
    //put every fit seq in own bin
    for (long long i = 0; i < randomizeSeqs.size(); i++) {
        // std::vector<long long> thisBin;
        // thisBin.emplace_back(randomizeSeqs[i]);
        bins.emplace_back(std::vector<long long>{randomizeSeqs[i]});
        seqBin[randomizeSeqs[i]] = numRefBins+i;

        const long long numCloseSeqs = (matrix->getNumFitClose(randomizeSeqs[i])); //does not include self
        fitfalseNegatives += numCloseSeqs;
        if (numCloseSeqs == 0) { numFitSingletons++; } //you are a singletons counted by the matrix as a fitSingleton, but you are not removed because you have ref dists we want to use in the fitting. Don't want to count you twice in stats output.
    }
    numFitSeqs = randomizeSeqs.size();

    fitfalseNegatives /= 2; //square matrix
    fittrueNegatives = numFitSeqs * (numFitSeqs-1)/2 - (fitfalsePositives + fitfalseNegatives + fittruePositives); //since everyone is a singleton no one clusters together. True negative = num far apart

    numComboSeqs = numRefSeqs + randomizeSeqs.size();

    combofalseNegatives = matrix->getNumDists() - reftruePositives; //number of distance in matrix for reference seqs - reftruePositives
    combotrueNegatives = numComboSeqs * (numComboSeqs-1)/2 - (reffalsePositives + reffalseNegatives + reftruePositives);
    combotruePositives = reftruePositives;
    combofalsePositives = reffalsePositives;

    const double comboValue = metric->getValue(combotruePositives, combotrueNegatives, combofalsePositives, combofalseNegatives);

    //add insert location
    seqBin[bins.size()] = -1;
    insertLocation = bins.size();
    // std::vector<long long> temp;
    bins.emplace_back();

    if (randomize) { Utils::mothurRandomShuffle(randomizeSeqs); }

    value = comboValue;
    // TODO: FInd why there are more slots in the bins than there are sequences...
    return value;

}
/***********************************************************************/
/* for each sequence with mutual information (close)
 * remove from current OTU and calculate MCC when sequence forms its own OTU or joins one of the other OTUs where there is a sequence within the `threshold` (no need to calculate MCC if the paired sequence is already in same OTU and no need to try every OTU - just those where there's a close sequence)
 * keep or move the sequence to the OTU where the `metric` is the largest - flip a coin on ties */
bool OptiFitCluster::update(double& listMetric) {
    //for each sequence (singletons removed on read)
    for (long long randomizeSeq : randomizeSeqs) {

        auto it = seqBin.find(randomizeSeq);

        int seqNumber = it->first;
        int binNumber = it->second;

        if (binNumber == -1) { }
        else {
            std::vector<long long> bestBin; bestBin.resize(2, binNumber);
            std::vector<double> tn; tn.push_back(fittrueNegatives); tn.push_back(combotrueNegatives);
            std::vector<double> tp; tp.push_back(fittruePositives); tp.push_back(combotruePositives);
            std::vector<double> fp; fp.push_back(fitfalsePositives); fp.push_back(combofalsePositives);
            std::vector<double> fn; fn.push_back(fitfalseNegatives); fn.push_back(combofalseNegatives);
            std::vector<double> bestMetric; bestMetric.resize(2, -1);  //bestMetric[0] = fitSeqs alone, bestMetric[1] = combo or ref and fit

            std::vector<double> bestTp; bestTp.resize(2, 0);
            std::vector<double> bestTn; bestTn.resize(2, 0);
            std::vector<double> bestFp; bestFp.resize(2, 0);
            std::vector<double> bestFn; bestFn.resize(2, 0);

            //close / far count in current bin
            std::vector<double> results = getCloseFarCounts(seqNumber, binNumber);
            double combocCount = results[0];  double combofCount = results[1];

            //close / far count in current bin for fit seqs
            std::vector<double> fitresults = getCloseFarFitCounts(seqNumber, binNumber);
            double fitcCount = fitresults[0];  double fitfCount = fitresults[1];

            //fit metrics in current bin
            bestMetric[0] = metric->getValue(tp[0], tn[0], fp[0], fn[0]);
            bestTp[0] = tp[0]; bestTn[0] = tn[0]; bestFp[0] = fp[0]; bestFn[0] = fn[0];

            //combo metric in current bin
            bestMetric[1] = metric->getValue(tp[1], tn[1], fp[1], fn[1]);
            bestTp[1] = tp[1]; bestTn[1] = tn[1]; bestFp[1] = fp[1]; bestFn[1] = fn[1];

            //if not already singleton, then calc value if singleton was created
            if (!((bins[binNumber].size()) == 1)) {
                //make a singleton
                fn[0]+=fitcCount; tn[0]+=fitfCount; fp[0]-=fitfCount; tp[0]-=fitcCount;
                fn[1]+=combocCount; tn[1]+=combofCount; fp[1]-=combofCount; tp[1]-=combocCount;

                double singleFitMetric = metric->getValue(tp[0], tn[0], fp[0], fn[0]);
                double singleComboMetric = metric->getValue(tp[1], tn[1], fp[1], fn[1]);
                if ((singleFitMetric > bestMetric[0]) || (singleComboMetric > bestMetric[1])) {
                    bestBin[1] = -1; bestTp[1] = tp[1]; bestTn[1] = tn[1]; bestFp[1] = fp[1]; bestFn[1] = fn[1];
                    bestMetric[1] = singleComboMetric;

                    bestBin[0] = -1; bestTp[0] = tp[0]; bestTn[0] = tn[0]; bestFp[0] = fp[0]; bestFn[0] = fn[0];
                    bestMetric[0] = singleFitMetric;
                }
            }

            std::set<long long> binsToTry;
            std::set<long long> closeSeqs = matrix->getCloseRefSeqs(seqNumber);
            for (long long closeSeq : closeSeqs) { binsToTry.insert(seqBin[closeSeq]); }

            //merge into each "close" otu
            std::vector<std::vector<double> > ties; std::vector<std::vector<double> > ties0;
            for (long long it : binsToTry) {
                //reset tn, tp,fp,fn values to original bin
                tn[0] = fittrueNegatives; tp[0] = fittruePositives; fp[0] = fitfalsePositives; fn[0] = fitfalseNegatives;
                tn[1] = combotrueNegatives; tp[1] = combotruePositives; fp[1] = combofalsePositives; fn[1] = combofalseNegatives;

                //move out of old bin
                fn[0]+=fitcCount; tn[0]+=fitfCount; fp[0]-=fitfCount; tp[0]-=fitcCount;
                fn[1]+=combocCount; tn[1]+=combofCount; fp[1]-=combofCount; tp[1]-=combocCount;

                results = getCloseFarCounts(seqNumber, it); //results[0] = close count, results[1] = far count
                fn[1]-=results[0]; tn[1]-=results[1];  tp[1]+=results[0]; fp[1]+=results[1]; //move into new bin

                results = getCloseFarFitCounts(seqNumber, it);
                fn[0]-=results[0]; tn[0]-=results[1];  tp[0]+=results[0]; fp[0]+=results[1]; //move into new bin - only consider fit seqs

                double newComboMetric = metric->getValue(tp[1], tn[1], fp[1], fn[1]); //score when sequence is moved
                double newFitMetric = metric->getValue(tp[0], tn[0], fp[0], fn[0]); //score when sequence is moved
                //new best
                if (newComboMetric > bestMetric[1]) {
                    ties.clear(); ties0.clear();
                    bestMetric[1] = newComboMetric; bestBin[1] = it; bestTp[1] = tp[1]; bestTn[1] = tn[1]; bestFp[1] = fp[1]; bestFn[1] = fn[1];
                    bestMetric[0] = newFitMetric; bestBin[0] = it; bestTp[0] = tp[0]; bestTn[0] = tn[0]; bestFp[0] = fp[0]; bestFn[0] = fn[0];
                    std::vector<double> tie; tie.push_back(bestMetric[1]); tie.push_back(bestBin[1]); tie.push_back(bestTp[1]);
                    tie.push_back(bestTn[1]); tie.push_back(bestFp[1]); tie.push_back(bestFn[1]); ties.push_back(tie);
                    std::vector<double> tie0; tie0.push_back(bestMetric[0]); tie0.push_back(bestBin[0]); tie0.push_back(bestTp[0]);
                    tie0.push_back(bestTn[0]); tie0.push_back(bestFp[0]); tie0.push_back(bestFn[0]); ties0.push_back(tie0);


                }else if (newComboMetric == bestMetric[1]) {
                    bestMetric[1] = newComboMetric; bestBin[1] = it; bestTp[1] = tp[1]; bestTn[1] = tn[1]; bestFp[1] = fp[1]; bestFn[1] = fn[1];
                    bestMetric[0] = newFitMetric; bestBin[0] = it; bestTp[0] = tp[0]; bestTn[0] = tn[0]; bestFp[0] = fp[0]; bestFn[0] = fn[0];
                    std::vector<double> tie; tie.push_back(bestMetric[1]); tie.push_back(bestBin[1]); tie.push_back(bestTp[1]);
                    tie.push_back(bestTn[1]); tie.push_back(bestFp[1]); tie.push_back(bestFn[1]); ties.push_back(tie);
                    std::vector<double> tie0; tie0.push_back(bestMetric[0]); tie0.push_back(bestBin[0]); tie0.push_back(bestTp[0]);
                    tie0.push_back(bestTn[0]); tie0.push_back(bestFp[0]); tie0.push_back(bestFn[0]); ties0.push_back(tie0);

                }
            }

            if (ties.size() > 1) {
                int randomTie = Utils::getRandomIndex((int)ties.size()-1);
                bestMetric[1] = ties[randomTie][0]; bestBin[1] = ties[randomTie][1]; bestTp[1] = ties[randomTie][2]; bestTn[1] = ties[randomTie][3]; bestFp[1] = ties[randomTie][4]; bestFn[1] = ties[randomTie][5];
                bestMetric[0] = ties0[randomTie][0]; bestBin[0] = ties0[randomTie][1]; bestTp[0] = ties0[randomTie][2]; bestTn[0] = ties0[randomTie][3]; bestFp[0] = ties0[randomTie][4]; bestFn[0] = ties0[randomTie][5];
            }

            //how to choose the best bin if they differ????
            long long newBin = bestBin[1];

            bool usedInsert = false;
            if (newBin == -1) {  newBin = insertLocation;  usedInsert = true;  }

            if (newBin != binNumber) {
                combotruePositives = bestTp[1]; combotrueNegatives = bestTn[1]; combofalsePositives = bestFp[1]; combofalseNegatives = bestFn[1];
                fittruePositives = bestTp[0]; fittrueNegatives = bestTn[0]; fitfalsePositives = bestFp[0]; fitfalseNegatives = bestFn[0];

                //move seq from i to j
                bins[newBin].push_back(seqNumber); //add seq to bestbin
                bins[binNumber].erase(remove(bins[binNumber].begin(), bins[binNumber].end(), seqNumber), bins[binNumber].end()); //remove from old bin i
            }

            if (usedInsert) { insertLocation = findInsert(); }

            //update seqBins
            seqBin[seqNumber] = newBin; //set new OTU location
        }
    }

    listMetric = metric->getValue(combotruePositives, combotrueNegatives, combofalsePositives, combofalseNegatives);
    return 0;

}
/***********************************************************************/
long long OptiFitCluster::getNumBins() {
    long long singletn = 0;

    singletn = matrix->getNumSingletons();

    for (const auto & bin : bins) { if (!bin.empty()) { singletn++; } }

    return singletn;
}
/***********************************************************************/
long long OptiFitCluster::getNumFitBins() {
    const ListVector* list = getFittedList("", false);

    int numBins = 0;
    if (list != nullptr) {
        numBins = list->getNumBins();
        delete list;
    }

    return numBins;
}

/***********************************************************************/
std::vector<double> OptiFitCluster::getStats(double& tp,  double& tn,  double& fp,  double& fn) {
    double singletn = 0;
    if (!closed) { singletn = matrix->getNumSingletons(); }
    double tempnumSeqs = numComboSeqs + singletn;

    tp = combotruePositives;
    fp = combofalsePositives;
    fn = combofalseNegatives;
    tn = tempnumSeqs * (tempnumSeqs-1)/2 - (combofalsePositives + combofalseNegatives + combotruePositives); //adds singletons to tn

    std::vector<double> results;

    Sensitivity sens;   double sensitivity = sens.getValue(tp, tn, fp, fn); results.push_back(sensitivity);
    Specificity spec;   double specificity = spec.getValue(tp, tn, fp, fn); results.push_back(specificity);
    PPV ppv;            double positivePredictiveValue = ppv.getValue(tp, tn, fp, fn); results.push_back(positivePredictiveValue);
    NPV npv;            double negativePredictiveValue = npv.getValue(tp, tn, fp, fn); results.push_back(negativePredictiveValue);
    FDR fdr;            double falseDiscoveryRate = fdr.getValue(tp, tn, fp, fn); results.push_back(falseDiscoveryRate);
    Accuracy acc;       double accuracy = acc.getValue(tp, tn, fp, fn); results.push_back(accuracy);
    MCC mcc;            double matthewsCorrCoef = mcc.getValue(tp, tn, fp, fn); results.push_back(matthewsCorrCoef);
    F1Score f1;         double f1Score = f1.getValue(tp, tn, fp, fn); results.push_back(f1Score);

    return results;
}
/***********************************************************************/
std::vector<double> OptiFitCluster::getFitStats(double& tp,  double& tn,  double& fp,  double& fn) {
    double singletn = 0;
    if (!closed) { singletn = matrix->getNumFitTrueSingletons(); }
    double tempnumSeqs = numFitSeqs + singletn; //numFitSingletons are reads that are selected as the fit seqs, that have dists to reference but not dists to other fit seqs. They are included

    tp = fittruePositives;
    fp = fitfalsePositives;
    fn = fitfalseNegatives;
    tn = tempnumSeqs * (tempnumSeqs-1)/2 - (fitfalsePositives + fitfalseNegatives + fittruePositives); //adds singletons to tn

    std::vector<double> results;

    Sensitivity sens;   double sensitivity = sens.getValue(tp, tn, fp, fn); results.push_back(sensitivity);
    Specificity spec;   double specificity = spec.getValue(tp, tn, fp, fn); results.push_back(specificity);
    PPV ppv;            double positivePredictiveValue = ppv.getValue(tp, tn, fp, fn); results.push_back(positivePredictiveValue);
    NPV npv;            double negativePredictiveValue = npv.getValue(tp, tn, fp, fn); results.push_back(negativePredictiveValue);
    FDR fdr;            double falseDiscoveryRate = fdr.getValue(tp, tn, fp, fn); results.push_back(falseDiscoveryRate);
    Accuracy acc;       double accuracy = acc.getValue(tp, tn, fp, fn); results.push_back(accuracy);
    MCC mcc;            double matthewsCorrCoef = mcc.getValue(tp, tn, fp, fn); results.push_back(matthewsCorrCoef);
    F1Score f1;         double f1Score = f1.getValue(tp, tn, fp, fn); results.push_back(f1Score);

    return results;
}
/***********************************************************************/
ListVector* OptiFitCluster::getList() {
    auto* list = new ListVector();
    const auto* singleton = matrix->getListSingle();

    if (singleton->size() <= 0) { //add in any sequences above cutoff in read. Removing these saves clustering time.
        for (int i = 0; i < singleton->getNumBins(); i++) {
            if (!singleton->get(i).empty()) {
                list->push_back(singleton->get(i));
            }
        }
        delete singleton;
    }

    for (const auto & bin : bins) {
        std::vector<int> thisBin;
        for (const long long j : bin) {  if (j >= 0) { thisBin.push_back(j); } }

        if (!thisBin.empty()) {

            std::string otu = matrix->getName(thisBin[0]);

            for (int j = 1; j < thisBin.size(); j++) { otu += "," + matrix->getName(thisBin[j]); }
            list->push_back(otu);
        }
    }

    return list;
}
/***********************************************************************/
ListVector* OptiFitCluster::getFittedList(std::string label, bool includerefs) {
    ListVector* list = new ListVector();

    std::unordered_map<long long, std::string> newBins;
    std::unordered_set<long long> unFitted;
    long long numListSeqs = 0;
    for (long long randomizeSeq : randomizeSeqs) { //build otus

        std::map<long long, long long>::iterator it = seqBin.find(randomizeSeq);

        long long seqNumber = it->first;
        long long binNumber = it->second;

        std::map<long long, std::string>::iterator itBinLabels = binLabels.find(binNumber); //do we have a label for this bin.  If the seq maps to existing bin then we should, otherwise we couldn't "fit" this sequence

        if (itBinLabels != binLabels.end()) {
            numListSeqs++;
            const auto itBin = newBins.find(binNumber); // have we seen this otu yet?

            if (itBin == newBins.end()) { //create bin
                newBins[binNumber] = matrix->getName(seqNumber);
            }else { //append bin
                newBins[binNumber] += "," + matrix->getName(seqNumber);
            }
        }else { unFitted.insert(seqNumber); }
    }

    if (denovo || includerefs) { //add in refs
        std::vector<long long> refs = matrix->getRefSeqs();

        for (long long ref : refs) {
            const auto it = seqBin.find(ref);

            long long seqNumber = it->first;
            long long binNumber = it->second;

            const auto itBin = newBins.find(binNumber); // have we seen this otu yet?

            if (itBin == newBins.end()) { //create bin
                newBins[binNumber] = matrix->getName(seqNumber);
            }else { //append bin
                newBins[binNumber] += "," + matrix->getName(seqNumber);
            }
        }
    }

    //numFitSeqs does not include any kind of singleton
    long long numUnFitted = (numFitSeqs + numFitSingletons - numListSeqs); //getNumFitTrueSingletons are fit reads that have no dists in the matrix. This can be confusing, think of it like this: there are true singletons, meaning we don't care if you are a ref or fit and you have no dists below the cutoff. This means you will be in your own OTU no matter what we do. There are fitSingletons, meaning you are a fit sequence and have no dists below the cutoff that coorespond to other fit seqs ( NOTE: you may or may not have dists to ref seqs or you could be a true singleton or a just a singleton because of the references chosen).

    long long numSingletonBins = 0;
    if ((label != "") && (numUnFitted != 0)) {

        // Rcpp::Rcout << ("\nFitted " + std::to_string(numListSeqs) + " sequences to " + std::to_string(newBins.size()) + " existing OTUs.\n");

        if (!closed) { //cluster the unfitted seqs separately
            // Rcpp::Rcout << (std::to_string(numUnFitted) + " sequences were unable to be fitted existing OTUs, excluding singletons.\n");

            // Rcpp::Rcout << ("\n**************** Clustering the unfitted sequences ****************\n");

            OptiData* unFittedMatrix = matrix->extractMatrixSubset(unFitted);
            //matrix->extractMatrixSubset(unFitted)

            ListVector* unfittedList = clusterUnfitted(unFittedMatrix, label); //unfittedList includes unfitted singletons

            if (unfittedList != nullptr) {

                // Rcpp::Rcout << ("The unfitted sequences clustered into " + std::to_string(unfittedList->getNumBins()) + " new OTUs.\n"); //+unFittedMatrix->getNumSingletons()+ matrix->getNumFitSingletons()

                for (int i = 0; i < unfittedList->getNumBins(); i++) {
                    std::string bin = unfittedList->get(i);
                    if (bin != "") { list->push_back(unfittedList->get(i)); }
                }
                delete unfittedList;
            }
            delete unFittedMatrix;

            // Rcpp::Rcout << ("\n*******************************************************************\n\n");

            //add in fit singletons
            ListVector* singleton = matrix->getFitListSingle();

            if (singleton != nullptr) { //add in any sequences above cutoff in read. Removing these saves clustering time.
                for (int i = 0; i < singleton->getNumBins(); i++) {
                    if (singleton->get(i) != "") { list->push_back(singleton->get(i)); }
                }
                numSingletonBins += singleton->getNumBins();
                delete singleton;
            }

        }else {
            // Rcpp::Rcout << ("\nSequences that were unable to be fitted existing OTUs will be listed in the *.optifit_scrap.accnos file.\n");
            unfittedNames = matrix->getNames(unFitted);

            //add in fit singletons
            ListVector* singleton = matrix->getFitListSingle();

            if (singleton != nullptr) { //add in any sequences above cutoff in read. Removing these saves clustering time.
                for (int i = 0; i < singleton->getNumBins(); i++) {
                    if (singleton->get(i) != "") { unfittedNames.insert(singleton->get(i)); }
                }
                delete singleton;
            }
        }
    }else {
        if (label != "") {
            // Rcpp::Rcout << ("\nFitted all " + std::to_string(list->getNumSeqs()) + " sequences to existing OTUs. \n");
        }
    }

    std::vector<std::string> newLabels = list->getLabels();

    for (auto &[fst, snd] : newBins) {
        list->push_back(snd);
        newLabels.push_back(binLabels[fst]);
    }

    list->setLabels(newLabels);

    return list;
}
/***********************************************************************/
int OptiFitCluster::findInsert() {

    //initially there are bins for each sequence (excluding singletons removed on read)
    for (int i = 0; i < bins.size(); i++) {
        if (bins[i].size() == 0) { return i;  } //this bin is empty
    }

    return -1;
}

/***********************************************************************/
ClusterExport* OptiFitCluster::runDenovoOptiCluster(std::map<std::string, int>& counts, std::string outStepFile){
    // Rcpp::message("\nClustering\n");
    // TODO: Reset optiflitcluster values
    OptifitClusterData *result = new OptifitClusterData("");
    const std::string cutoffString = std::to_string(cutoff);
    OptifitClusterInformation clusterInformation;
    clusterInformation.numberOfOtu = -1;
    constexpr double stableMetric = 0;
    constexpr int maxIters = 100;
    bool printStepsHeader = true;
    constexpr int denovoIters = 100;
    int smallestBins = -1;
    size_t index = 0;
    // Could thread?
    for (int i = 0; i < denovoIters; i++) {

        // OptiFitCluster cluster(matrix, metric, 0);
        // tag = cluster.getTag();
        int iters = 0;
        double listVectorMetric = 0; //worst state
        double delta = 1;

        //get "ref" seqs for initialize inputs
        OptiData* refMatrix = matrix->extractRefMatrix(); // This is an empty matrix, that is the issue
        OptiCluster clust(refMatrix, new MCC(), cutoff, 0.0001, 0);
        const auto exportResult = clust.Execute();//clusterRefs(refMatrix, metric);
        ListVector refList = exportResult->GetListVector().listVector;
        delete refMatrix;

        std::vector<std::vector<std::string> > otus;
        for (int j = 0; j < refList.getNumBins(); j++) {
            if (const std::string bin = refList.get(j); !bin.empty()) {
                std::vector<std::string> binNames;
                Utils::splitAtComma(bin, binNames);
                otus.push_back(binNames);
            }
        }

        //add tag to OTULabels to indicate the reference
        std::vector<std::string> refListLabels = refList.getLabels();
        for (auto & refListLabel : refListLabels) { refListLabel = "Ref_" + refListLabel;  }
        refList.setLabels(refListLabels);

        initialize(listVectorMetric, true, otus, refList.getLabels(), "", true);

        // delete refList;

        long long numBins = getNumBins();

        double tp, tn, fp, fn;
        std::string clusterMetrics;
        std::string sensFile;
        std::vector<double> stats = getStats(tp, tn, fp, fn);
        std::vector<std::string> clusterMetricList;
        double fittp, fittn, fitfp, fitfn;
        long long numFitBins = getNumFitBins();
        std::vector<double> fitresults = getFitStats(fittp, fittn, fitfp, fitfn);

        // Rcpp::message("\nFitting " + std::to_string(matrix->getNumFitSeqs()+matrix->getNumFitSingletons()+matrix->getNumFitTrueSingletons()) + " sequences to reference otus.\n");

        // Rcpp::message("\n\nlist\tstate\titer\tlabel\tnum_otus\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n");

        std::vector<std::string> sensfileHeaders{"label","cutoff","tp","tn","fp","fn","sensitivity",
            "specificity","ppv","npv","fdr","accuracy","mcc","f1score"};
        // sensFile = "label\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n";

        std::vector<std::string> clusterMetricsHeaders{"iter","time","label","num_otus","cutoff","tp","tn",
            "fp","fn","sensitivity","specificity","ppv","npv",
            "fdr", "accuracy", "mcc", "f1score"};

        // clusterMetrics = ("0,0," + std::to_string(cutoff) + "," + std::to_string(numBins) + "," +
        //                       std::to_string(cutoff) + "," + std::to_string(tp) + "," + std::to_string(tn) + "," +
        //                       std::to_string(fp) + "," + std::to_string(fn) + ",");
        // for (double stat: stats) {
        //     clusterMetrics += (std::to_string(stat) + ",");
        // }

        // outputSteps(outStepFile, printStepsHeader, tp, tn, fp, fn, results, numBins, fittp, fittn, fitfp, fitfn, fitresults, numFitBins, 0, false, 0);

        while ((delta > stableMetric) && (iters < maxIters)) { //

            double oldMetric = listVectorMetric;
            auto startTime = std::chrono::system_clock::now();
            update(listVectorMetric);

            delta = abs(oldMetric - listVectorMetric);
            iters++;

            stats = getStats(tp, tn, fp, fn);
            ListVector* listVector = getFittedList(cutoffString, false);
            numFitBins = listVector->getNumBins();
            delete listVector;
            fitresults = getFitStats(fittp, fittn, fitfp, fitfn);
            numBins = getNumBins();
            auto endTime = std::chrono::system_clock::now();
            std::chrono::duration<double> currentTime = endTime - startTime;
            clusterMetrics = (std::to_string(i) + "," + std::to_string(currentTime.count()) + "," +
                              cutoffString + "," + std::to_string(numFitBins) + "," +
                              cutoffString + "," + std::to_string(fittp) + "," + std::to_string(fittn) + ","
                              + std::to_string(fitfp) + "," + std::to_string(fitfn) + ",");
            for (double stat: fitresults) {
                clusterMetrics += (std::to_string(stat) + ",");
            }
            Utils::AddRowToDataFrameMap(dataframeMapClusterMetrics, clusterMetrics, clusterMetricsHeaders);

            // outputSteps(outStepFile, printStepsHeader, tp, tn, fp, fn, results, numBins, fittp, fittn, fitfp, fitfn, fitresults, numFitBins, iters, false, i);
        }
        // outputSteps(outStepFile, printStepsHeader, tp, tn, fp, fn, results, numBins, fittp, fittn, fitfp, fitfn, fitresults, numFitBins, iters, true, i);
        //m->mothurOutEndLine(); m->mothurOutEndLine();
        // f (m->getControl_pressed()) {  return 0; }
        //
        // ofstream listFile;
        // tag = "optifit_" + metric->getName() + "_denovo." + toString(i+1);
        // string listFileName = fileroot+ tag + ".list";
        // util.openOutputFile(listFileName,    listFile);

        ListVector* list = getFittedList(cutoffString, false);
        // if (clusterInformation.numberOfOtu < list->getNumBins()) {
        //     delete list;
        //     continue;
        // }

        if (const int binCount = static_cast<int>(list->getNumBins());
            clusterInformation.numberOfOtu == -1 || binCount < clusterInformation.numberOfOtu) {
            clusterInformation.label = cutoffString;
            clusterInformation.numberOfOtu = binCount;
            // clusterInformation.clusterBins = list->print(listFile);
            // Rcpp::Rcout << clusterInformation.clusterBins << std::endl;
            result->SetListVector(*list, cutoffString);
            // result->GetListVector().listVector.print()
            list->setLabel(std::to_string(cutoff));
        }

        sensFile += cutoffString + ',' + cutoffString + ',' + std::to_string(tp) + ',' +
                std::to_string(tn) + ',' +
                std::to_string(fp) + ',' + std::to_string(fn) + ',';
        for (double res: stats) { sensFile += std::to_string(res) + ','; }
        Utils::AddRowToDataFrameMap(dataframeMapSensMetrics, sensFile, sensfileHeaders);
        // list->setLabels(nullVector);

        // if(countfile != "") { list->print(listFile, counts); }
        // else { list->print(listFile); }

        // listFile.close();
        // listFiles.push_back(listFileName);
        // if (smallestBins == -1 || smallestBins > binCount) {
        //     smallestBins = binCount;
        //     index = i;
        // }

        delete list;
        Reset();
        matrix->randomizeRefs();
    }

    // tag = "optifit_" + metric->getName() + "_denovo";
    // string listFileName = fileroot+ tag + ".list";
    std::ofstream listFile;
    clusterInformation.clusterBins = result->GetListVector().listVector.print(listFile);
    result->AddToData(clusterInformation);
    return result;
}
/***********************************************************************/
std::vector<double> OptiFitCluster::getCloseFarCounts(const long long seq, const long long newBin) {
    std::vector<double> results; results.push_back(0); results.push_back(0); //results[0] = close count, results[1] = far count

    if (newBin == -1) { }  //making a singleton bin. Close but we are forcing apart.
    else { //merging a bin
        for (const long long currentBinSeq : bins[newBin]) {
            if (seq == currentBinSeq) {} //ignore self
            else if (!matrix->isClose(seq, currentBinSeq)) { results[1]++; }  //this sequence is "far away" from sequence i - above the cutoff
            else { results[0]++;  }  //this sequence is "close" to sequence i - distance between them is less than cutoff
        }
    }

    return results;
}
/***********************************************************************/
std::vector<double> OptiFitCluster::getCloseFarFitCounts(long long seq, long long newBin) {
    std::vector<double> results; results.push_back(0); results.push_back(0); //results[0] = close count, results[1] = far count

    if (newBin == -1) { }  //making a singleton bin. Close but we are forcing apart.
    else { //merging a bin
        for (const long long currentBinSeq : bins[newBin]) {

            if (seq == currentBinSeq) {} //ignore self
            else {
                bool isFit = true;
                if (bool closeFit = matrix->isCloseFit(seq, currentBinSeq, isFit)) { //you are close if you are fit and close
                    results[0]++;
                }else if (isFit) { results[1]++; } //this sequence is "far away" and fit - above the cutoff
            }
        }
    }
    return results;
}

/***********************************************************************/
ListVector* OptiFitCluster::clusterUnfitted(OptiData* unfittedMatrix, std::string label) {

    ListVector* list = nullptr;
    auto* opti = new OptiMatrix();
    opti->MoveData(unfittedMatrix);
    OptiCluster cluster(opti, metric, opti->GetCutoff(), 0);
    int iters = 0;
    double listVectorMetric = 0; //worst state
    double delta = 1;

    cluster.initialize(listVectorMetric, true, "singleton");

    long long numBins = cluster.getNumBins();
    // Rcpp::Rcout << ("\n\niter\ttime\tlabel\tnum_otus\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n");

    double tp, tn, fp, fn;
    std::vector<double> results = cluster.getStats(tp, tn, fp, fn);
    // Rcpp::Rcout << ("0\t0\t" + label + "\t" + std::to_string(numBins) + "\t"+ label + "\t" + std::to_string(tp) + "\t" + std::to_string(tn) + "\t" + std::to_string(fp) + "\t" + std::to_string(fn) + "\t");

    // for (const double result : results) { Rcpp::Rcout << (std::to_string(result) + "\t");  }

    while ((delta > 0.0001) && (iters < 100)) {

        const long start = time(nullptr);

        const double oldMetric = listVectorMetric;

        cluster.update(listVectorMetric);

        delta = abs(oldMetric - listVectorMetric);
        iters++;

        results = cluster.getStats(tp, tn, fp, fn);
        numBins = cluster.getNumBins();

        // Rcpp::Rcout << (std::to_string(iters) + "\t" + std::to_string(time(nullptr) - start) + "\t" + label + "\t" + std::to_string(numBins) + "\t" + label + "\t"+ std::to_string(tp) + "\t" + std::to_string(tn) + "\t" + std::to_string(fp) + "\t" + std::to_string(fn) + "\t");

        // for (const double result : results) { Rcpp::Rcout << (std::to_string(result) + "\t");  }

    }
    list = new ListVector(cluster.getList());
    list->setLabel(label);

    return list;
}

// ClusterExport OptiFitCluster::runUserRefOptiCluster(OptiData*& matrix, ClusterMetric*& metric, map<string, int>& counts, string outStepFile, vector<string> refListLabels, vector<vector<string> > otus){
//
//         bool printStepsHeader = true;
//         // tag = cluster.getTag();
//
//         int iters = 0;
//         double listVectorMetric = 0; //worst state
//         double delta = 1;
//
//         if (!createAccnos) {
//
//             m->mothurOut("\nClustering references from " + distfile + "\n");
//
//             //get "ref" seqs for initialize inputs
//             OptiData* refMatrix = matrix->extractRefMatrix();
//             ListVector* refList = clusterRefs(refMatrix, metric); delete refMatrix;
//
//             for (int i = 0; i < refList->getNumBins(); i++) {
//                 vector<string> binNames;
//                 string bin = refList->get(i);
//                 if (bin != "") {
//                     util.splitAtComma(bin, binNames);
//                     otus.push_back(binNames);
//                 }
//             }
//
//             //add tag to OTULabels to indicate the reference
//             refListLabels = refList->getLabels();
//             for (int i = 0; i < refListLabels.size(); i++) { refListLabels[i] = "Ref_" + refListLabels[i];  }
//             refList->setLabels(refListLabels);
//             delete refList;
//         }
//
//         cluster.initialize(listVectorMetric, true, otus, refListLabels, method, false);
//
//         long long numBins = cluster.getNumBins();
//         double tp, tn, fp, fn;
//         vector<double> results = cluster.getStats(tp, tn, fp, fn);
//
//         double fittp, fittn, fitfp, fitfn;
//         long long numFitBins = cluster.getNumFitBins();
//         vector<double> fitresults = cluster.getFitStats(fittp, fittn, fitfp, fitfn);
//
//         m->mothurOut("\nFitting " + toString(matrix->getNumFitSeqs()+matrix->getNumFitSingletons()+matrix->getNumFitTrueSingletons()) + " sequences to reference otus.\n");
//
//         m->mothurOut("\n\nlist\tstate\titer\tlabel\tnum_otus\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n");
//
//         outputSteps(outStepFile, printStepsHeader, tp, tn, fp, fn, results, numBins, fittp, fittn, fitfp, fitfn, fitresults, numFitBins, 0, false, 0);
//
//
//         while ((delta > stableMetric) && (iters < maxIters)) { //
//
//             if (m->getControl_pressed()) { break; }
//             double oldMetric = listVectorMetric;
//
//             cluster.update(listVectorMetric);
//
//             delta = abs(oldMetric - listVectorMetric);
//             iters++;
//
//             results = cluster.getStats(tp, tn, fp, fn);
//             numBins = cluster.getNumBins();
//             numFitBins = cluster.getNumFitBins();
//             fitresults = cluster.getFitStats(fittp, fittn, fitfp, fitfn);
//
//             outputSteps(outStepFile, printStepsHeader, tp, tn, fp, fn, results, numBins, fittp, fittn, fitfp, fitfn, fitresults, numFitBins, iters, false, 0);
//         }
//         m->mothurOutEndLine(); m->mothurOutEndLine();
//
//         if (m->getControl_pressed()) {  return 0; }
//
//         ListVector* list = cluster.getFittedList(toString(cutoff), printref);
//         list->setLabel(toString(cutoff));
//
//         string sensspecFilename = fileroot+ tag + ".sensspec";
//         ofstream sensFile;
//         util.openOutputFile(sensspecFilename,    sensFile);
//         outputNames.push_back(sensspecFilename); outputTypes["sensspec"].push_back(sensspecFilename);
//
//         if (method == "closed") {
//             sensFile << "label\tcutoff\tnumotus\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n";
//             int numBins = list->getNumBins();
//             if (printref) { //combo
//                 results = cluster.getStats(tp, tn, fp, fn);
//                 sensFile << cutoff << '\t' << cutoff << '\t' << numBins << '\t' << tp << '\t' << tn << '\t' << fp << '\t' << fn;
//                 for (int i = 0; i < results.size(); i++) {  sensFile << '\t' << results[i]; } sensFile << '\n';
//             }else { //fit
//                 fitresults = cluster.getFitStats(fittp, fittn, fitfp, fitfn);
//                 sensFile << cutoff << '\t' << cutoff << '\t' << numBins << '\t' << fittp << '\t' << fittn << '\t' << fitfp << '\t' << fitfn;
//                 for (int i = 0; i < fitresults.size(); i++) {  sensFile << "\t" << fitresults[i]; } sensFile << endl;
//             }
//             set<string> unfitted = cluster.getUnfittedNames();
//
//             string accnosFilename = fileroot+ "optifit_scrap.accnos";
//             outputNames.push_back(accnosFilename); outputTypes["accnos"].push_back(accnosFilename);
//
//             ofstream accOut; util.openOutputFile(accnosFilename,    accOut);
//             for (set<string>::iterator it = unfitted.begin(); it != unfitted.end(); it++) {
//                 accOut << *it << endl;
//             }
//             accOut.close();
//         }else {
//             runSensSpec(matrix, metric, list, counts, sensFile);
//         }
//         sensFile.close();
//
//         return list;
// }



ClusterExport* OptiFitCluster::runUserRefOptiCluster(ClusterMetric*& metric, std::vector<std::string> refListLabels,
std::vector<std::vector<std::string> > otus){
    bool printStepsHeader = true;
    const std::string cutoffString = std::to_string(cutoff);
    constexpr double stableMetric = 0;
    constexpr int maxIters = 100;
    OptifitClusterData *result = new OptifitClusterData("");
    OptifitClusterInformation clusterInformation;

    // OptiFitCluster cluster(matrix, metric, cutoff, 0);
    // tag = cluster.getTag();

    int iters = 0;
    double listVectorMetric = 0; //worst state
    double delta = 1;

// if (!createAccnos) {

    //m->mothurOut("\nClustering references from " + distfile + "\n");

    //get "ref" seqs for initialize inputs
    // OptiData* refMatrix = matrix->extractRefMatrix();
    // ListVector refList = clusterRefs(refMatrix, metric); delete refMatrix;
    //
    // for (int i = 0; i < refList.getNumBins(); i++) {
    //     std::vector<std::string> binNames;
    //     std::string bin = refList.get(i);
    //     if (bin != "") {
    //         Utils::splitAtComma(bin, binNames);
    //         otus.push_back(binNames);
    //     }
    // }
    //
    // // add tag to OTULabels to indicate the reference
    //  refListLabels = refList.getLabels();
    //  for (auto & refListLabel : refListLabels) { refListLabel = "Ref_" + refListLabel;  }
    //  refList.setLabels(refListLabels);
 //}

    initialize(listVectorMetric, true, otus, refListLabels, "", false);

    long long numBins = getNumBins();
    double tp, tn, fp, fn;
    std::vector<double> results = getStats(tp, tn, fp, fn);

    double fittp, fittn, fitfp, fitfn;
    long long numFitBins = getNumFitBins();
    std::vector<double> fitresults; //= getFitStats(fittp, fittn, fitfp, fitfn);
    std::string clusterMetrics;
    std::string sensFile;
    // std::vector<double> stats = getStats(tp, tn, fp, fn);
    std::vector<std::string> clusterMetricList;

    // Rcpp::message("\nFitting " + std::to_string(matrix->getNumFitSeqs()+matrix->getNumFitSingletons()+matrix->getNumFitTrueSingletons()) + " sequences to reference otus.\n");

    // Rcpp::message("\n\nlist\tstate\titer\tlabel\tnum_otus\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n");

    std::vector<std::string> sensfileHeaders{"label","cutoff","tp","tn","fp","fn","sensitivity",
    "specificity","ppv","npv","fdr","accuracy","mcc","f1score"};
    // sensFile = "label\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n";

    std::vector<std::string> clusterMetricsHeaders{"iter","time","label","num_otus","cutoff","tp","tn",
      "fp","fn","sensitivity","specificity","ppv","npv",
        "fdr", "accuracy", "mcc", "f1score"};
    // m->mothurOut("\nFitting " + std::to_string(matrix->getNumFitSeqs()+matrix->getNumFitSingletons()+matrix->getNumFitTrueSingletons()) + " sequences to reference otus.\n");

    // Rcpp::message("\n\nlist\tstate\titer\tlabel\tnum_otus\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n");

    // outputSteps(outStepFile, printStepsHeader, tp, tn, fp, fn, results, numBins, fittp, fittn, fitfp, fitfn, fitresults, numFitBins, 0, false, 0);


    while ((delta > stableMetric) && (iters < maxIters)) { //

        // if (m->getControl_pressed()) { break; }
        const double oldMetric = listVectorMetric;
        auto startTime = std::chrono::system_clock::now();
        update(listVectorMetric);

        delta = abs(oldMetric - listVectorMetric);
        iters++;

        results = getStats(tp, tn, fp, fn);
        numBins = getNumBins();
        numFitBins = getNumFitBins();
        fitresults = getFitStats(fittp, fittn, fitfp, fitfn);
        auto endTime = std::chrono::system_clock::now();
        std::chrono::duration<double> currentTime = endTime - startTime;
        clusterMetrics = (std::to_string(iters) + "," + std::to_string(currentTime.count()) + "," +
                           cutoffString + "," + std::to_string(numFitBins) + "," +
                           cutoffString + "," + std::to_string(fittp) + "," + std::to_string(fittn) + ","
                           + std::to_string(fitfp) + "," + std::to_string(fitfn) + ",");
        for (double stat: fitresults) {
            clusterMetrics += (std::to_string(stat) + ",");
        }
        Utils::AddRowToDataFrameMap(dataframeMapClusterMetrics, clusterMetrics, clusterMetricsHeaders);

        // outputSteps(outStepFile, printStepsHeader, tp, tn, fp, fn, results, numBins, fittp, fittn, fitfp, fitfn, fitresults, numFitBins, iters, false, 0);
    }
    //    m->mothurOutEndLine(); m->mothurOutEndLine();

    // if (m->getControl_pressed()) {  return 0; }
    std::ofstream listFile;
    ListVector* list = getFittedList(std::to_string(cutoff), true);
    list->setLabel(std::to_string(cutoff));
    clusterInformation.label = std::to_string(cutoff);
    clusterInformation.numberOfOtu = static_cast<int>(numBins);
    clusterInformation.clusterBins = list->print(listFile);
    result->AddToData(clusterInformation);
    result->SetListVector(*list, std::to_string(cutoff));
    delete list;

    // string sensspecFilename = fileroot+ tag + ".sensspec";
    // ofstream sensFile;
    // util.openOutputFile(sensspecFilename,    sensFile);
    // outputNames.push_back(sensspecFilename); outputTypes["sensspec"].push_back(sensspecFilename);

    // if (method == "closed") {
    //     sensFile << "label\tcutoff\tnumotus\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n";
    //     int numBins = list->getNumBins();
    //     if (printref) { //combo
    //         results = cluster.getStats(tp, tn, fp, fn);
    //         sensFile << cutoff << '\t' << cutoff << '\t' << numBins << '\t' << tp << '\t' << tn << '\t' << fp << '\t' << fn;
    //         for (int i = 0; i < results.size(); i++) {  sensFile << '\t' << results[i]; } sensFile << '\n';
    //     }else { //fit
    //         fitresults = cluster.getFitStats(fittp, fittn, fitfp, fitfn);
    //         sensFile << cutoff << '\t' << cutoff << '\t' << numBins << '\t' << fittp << '\t' << fittn << '\t' << fitfp << '\t' << fitfn;
    //         for (int i = 0; i < fitresults.size(); i++) {  sensFile << "\t" << fitresults[i]; } sensFile << endl;
    //     }
    //     set<string> unfitted = cluster.getUnfittedNames();
    //
    //     string accnosFilename = fileroot+ "optifit_scrap.accnos";
    //     outputNames.push_back(accnosFilename); outputTypes["accnos"].push_back(accnosFilename);
    //
    //     ofstream accOut; util.openOutputFile(accnosFilename,    accOut);
    //     for (set<string>::iterator it = unfitted.begin(); it != unfitted.end(); it++) {
    //         accOut << *it << endl;
    //     }
    //     accOut.close();
    // }else {
    //     runSensSpec(matrix, metric, list, counts, sensFile);
    // }
    // sensFile.close();

    return result;
}



ListVector OptiFitCluster::clusterRefs(OptiData*& refsMatrix, ClusterMetric*& metric) {

    // Rcpp::message("\nClustering " + std::to_string(refsMatrix->GetNameList().size()+refsMatrix->getNumSingletons()) + " reference sequences.\n");

    ListVector list;

    OptiCluster cluster(refsMatrix, metric, 0, 0);
    constexpr int maxIters = 100;
    int iters = 0;
    double listVectorMetric = 0; //worst state
    double delta = 1;

    cluster.initialize(listVectorMetric, true, "singleton");

    long long numBins = cluster.getNumBins();
    // Rcpp::message("\n\niter\ttime\tlabel\tnum_otus\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n");

    double tp, tn, fp, fn;
    std::vector<double> results = cluster.getStats(tp, tn, fp, fn);
    // Rcpp::message("0\t0\t" + std::to_string(cutoff) + "\t" + std::to_string(numBins) + "\t"+ std::to_string(cutoff) +
    //     "\t" + std::to_string(tp) + "\t" + std::to_string(tn) + "\t" + std::to_string(fp) + "\t" + std::to_string(fn) +
    //     "\t");

    for (const double result : results) {
        // Rcpp::message(std::to_string(result) + "\t");
    }

    while ((delta > 0.0001) && (iters < maxIters)) {

        long start = time(nullptr);

        double oldMetric = listVectorMetric;

        cluster.update(listVectorMetric);

        delta = abs(oldMetric - listVectorMetric);
        iters++;

        results = cluster.getStats(tp, tn, fp, fn);
        numBins = cluster.getNumBins();

        // Rcpp::message(std::to_string(iters) + "\t" + std::to_string(time(nullptr) - start) + "\t" +
        //     std::to_string(cutoff) + "\t" + std::to_string(numBins) + "\t" + std::to_string(cutoff) + "\t"+
        //     std::to_string(tp) + "\t" + std::to_string(tn) + "\t" + std::to_string(fp) + "\t" + std::to_string(fn) +
        //     "\t");

        // for (const double result : results) { Rcpp::message(std::to_string(result) + "\t");  }

    }
    list = cluster.getList();
    list.setLabel(std::to_string(cutoff));

    return list;
}

ClusterExport* OptiFitCluster::runRefOptiCluster(ListVector refList, std::map<std::string, int>& counts, std::string outStepFile){
    int iters = 0;
    double listVectorMetric = 0; //worst state
    double delta = 1;
    int maxIters = 100;
    constexpr double stableMetric = 0;

    std::vector<std::vector<std::string> > otus;
    for (int i = 0; i < refList.getNumBins(); i++) {
        std::vector<std::string> binNames;
        if (const std::string bin = refList.get(i); !bin.empty()) {
            Utils::splitAtComma(bin, binNames);
            otus.push_back(binNames);
        }
    }

    std::map<std::string, int> refCounts;
    // if (refcountfile != "") {
    //     CountTable refct; refct.readTable(refcountfile, false, false);
    //     refCounts = refct.getNameMap();
    // }else if (refnamefile != "") { refCounts = util.readNames(refnamefile); }
    // else { //assume unique
    //     for (auto & otu : otus) { for (int j = 0; j < otu.size(); j++) { refCounts[otu[j]] = 1; } }
    // }
    counts.insert(refCounts.begin(), refCounts.end());

    initialize(listVectorMetric, true, otus, refList.getLabels(), "", false);
    const std::string cutoffString = std::to_string(cutoff);
    long long numBins = getNumBins();
    double tp, tn, fp, fn;
    std::vector<double> results = getStats(tp, tn, fp, fn);
    std::string clusterMetrics;
    std::string sensFile;
    std::vector<double> stats = getStats(tp, tn, fp, fn);
    std::vector<std::string> clusterMetricList;
    double fittp, fittn, fitfp, fitfn;
    long long numFitBins = getNumFitBins();
    std::vector<double> fitresults = getFitStats(fittp, fittn, fitfp, fitfn);

    bool printStepsHeader = true;
    // outputSteps(outStepFile, printStepsHeader, tp, tn, fp, fn, results, numBins, fittp, fittn, fitfp, fitfn, fitresults, numFitBins, 0, true, 0);

    std::vector<std::string> sensfileHeaders{"label","cutoff","tp","tn","fp","fn","sensitivity",
  "specificity","ppv","npv","fdr","accuracy","mcc","f1score"};
    // sensFile = "label\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n";

    std::vector<std::string> clusterMetricsHeaders{"iter","time","label","num_otus","cutoff","tp","tn",
      "fp","fn","sensitivity","specificity","ppv","npv",
        "fdr", "accuracy", "mcc", "f1score"};
    while ((delta > stableMetric) && (iters < maxIters)) { //

       //  if (m->getControl_pressed()) { break; }
        auto startTime = std::chrono::system_clock::now();
        double oldMetric = listVectorMetric;

        update(listVectorMetric);

        delta = abs(oldMetric - listVectorMetric);
        iters++;

        // results = getStats(tp, tn, fp, fn);
        numBins = getNumBins();
        numFitBins = getNumFitBins();
        fitresults = getFitStats(fittp, fittn, fitfp, fitfn);
        auto endTime = std::chrono::system_clock::now();
        std::chrono::duration<double> currentTime = endTime - startTime;
        clusterMetrics = (std::to_string(iters) + "," + std::to_string(currentTime.count()) + "," +
                           cutoffString + "," + std::to_string(numFitBins) + "," +
                           cutoffString + "," + std::to_string(fittp) + "," + std::to_string(fittn) + ","
                           + std::to_string(fitfp) + "," + std::to_string(fitfn) + ",");
        for (double stat: fitresults) {
            clusterMetrics += (std::to_string(stat) + ",");
        }
        Utils::AddRowToDataFrameMap(dataframeMapClusterMetrics, clusterMetrics, clusterMetricsHeaders);


        // outputSteps(outStepFile, printStepsHeader, tp, tn, fp, fn, results, numBins, fittp, fittn, fitfp, fitfn, fitresults, numFitBins, iters, true, 0);
    }
    // m->mothurOutEndLine(); m->mothurOutEndLine();

    // if (m->getControl_pressed()) {  return 0; }
    std::ofstream listFile;
    auto* result = new OptifitClusterData("");
    ListVector* list = getFittedList(std::to_string(cutoff), true);
    list->setLabel(std::to_string(cutoff));
    OptifitClusterInformation clusterInformation;
    clusterInformation.label = std::to_string(cutoff);
    clusterInformation.numberOfOtu = static_cast<int>(numBins);
    clusterInformation.clusterBins = list->print(listFile);
    result->AddToData(clusterInformation);
    result->SetListVector(*list, std::to_string(cutoff));
    list->setLabel(std::to_string(cutoff));

    // ofstream listFile;
    // string listFileName = fileroot+ tag + ".list";
    // util.openOutputFile(listFileName,    listFile);

    // if(countfile != "") { list->print(listFile, counts); }
    // else { list->print(listFile); }
    // listFile.close();
    //
    // string sensspecFilename = fileroot+ tag + ".sensspec";
    // ofstream sensFile;
    // util.openOutputFile(sensspecFilename,    sensFile);
    // outputNames.push_back(sensspecFilename); outputTypes["sensspec"].push_back(sensspecFilename);
    //
    // if (method == "closed") {
    //     sensFile << "label\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n";
    //
    //     if (printref) { //combo
    //         results = cluster.getStats(tp, tn, fp, fn);
    //         sensFile << cutoff << '\t' << cutoff << '\t' << tp << '\t' << tn << '\t' << fp << '\t' << fn;
    //         for (int i = 0; i < results.size(); i++) {  sensFile << '\t' << results[i]; } sensFile << '\n';
    //     }else { //fit
    //         fitresults = cluster.getFitStats(fittp, fittn, fitfp, fitfn);
    //         sensFile << cutoff << '\t' << cutoff << '\t' << fittp << '\t' << fittn << '\t' << fitfp << '\t' << fitfn;
    //         for (int i = 0; i < fitresults.size(); i++) {  sensFile << "\t" << fitresults[i]; } sensFile << endl;
    //     }
    //     set<string> unfitted = cluster.getUnfittedNames();
    //
    //     string accnosFilename = fileroot+ "optifit_scrap.accnos";
    //     outputNames.push_back(accnosFilename); outputTypes["accnos"].push_back(accnosFilename);
    //
    //     ofstream accOut; util.openOutputFile(accnosFilename,    accOut);
    //     for (set<string>::iterator it = unfitted.begin(); it != unfitted.end(); it++) {
    //         accOut << *it << endl;
    //     }
    //     accOut.close();
    //
    // }else {
    //     runSensSpec(matrix, metric, list, counts, sensFile);
    // }
   // sensFile.close();

    delete list;
    return result;
    // return listFileName;

}