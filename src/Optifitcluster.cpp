

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
#include "DataStructures/OptiRefMatrix.h"


/***********************************************************************/
OptiFitCluster::OptiFitCluster(OptiRefMatrix* mt, ClusterMetric* met, const long long ns)
 : metric(met), matrix(mt), numComboSingletons(ns) {
    maxRefBinNumber = 0;
    closed = false;

    numFitSeqs = 0;  fittruePositives = 0; fitfalsePositives = 0; fitfalseNegatives = 0; fittrueNegatives = 0; numFitSingletons = 0;
    numComboSeqs = 0; numComboSingletons = 0; combotruePositives = 0; combofalsePositives = 0; combofalseNegatives = 0; combotrueNegatives = 0;
}

ClusterExport * OptiFitCluster::Execute() {
       if (abort) { if (calledHelp) { return 0; }  return 2;    }

        time_t estart = time(nullptr);
        
        ClusterMetric* metric = nullptr;
        if (metricName == "mcc")             { metric = new MCC();              }
        else if (metricName == "sens")       { metric = new Sensitivity();      }
        else if (metricName == "spec")       { metric = new Specificity();      }
        else if (metricName == "tptn")       { metric = new TPTN();             }
        else if (metricName == "tp")         { metric = new TP();               }
        else if (metricName == "tn")         { metric = new TN();               }
        else if (metricName == "fp")         { metric = new FP();               }
        else if (metricName == "fn")         { metric = new FN();               }
        else if (metricName == "f1score")    { metric = new F1Score();          }
        else if (metricName == "accuracy")   { metric = new Accuracy();         }
        else if (metricName == "ppv")        { metric = new PPV();              }
        else if (metricName == "npv")        { metric = new NPV();              }
        else if (metricName == "fdr")        { metric = new FDR();              }
        else if (metricName == "fpfn")       { metric = new FPFN();             }

        map<string, int> counts;
        string dupsFile = countfile; nameOrCount = "count";
        if (namefile != "") { dupsFile = namefile; nameOrCount = "name"; }
        else { CountTable ct; ct.readTable(countfile, false, false); counts = ct.getNameMap();  }

        if (outputdir == "") { outputdir += util.hasPath(distfile); }
        fileroot = outputdir + util.getRootName(util.getSimpleName(distfile));

        string listFile = ""; string bestListFileName = ""; string outputName = "";

        if (selfReference) { //de novo

            map<string, string> variables;
            variables["[filename]"] = fileroot;
            variables["[clustertag]"] = "optifit_" + metric->getName();
            outputName = getOutputFileName("steps", variables);

            if ((accnosfile == "") && (!createAccnos)) { //denovo with mothur randomly assigning references

                m->mothurOut("\nRandomly assigning reads from " + distfile + " as reference sequences\n");

                //distfile, distFormat, dupsFile, dupsFormat, cutoff, percentage to be fitseqs - will randomly assign as fit
                OptiData* matrix = new OptiRefMatrix(distfile, "column", dupsFile, nameOrCount, cutoff, fitPercent, refWeight);

                runDenovoOptiCluster(matrix, metric, counts, outputName);

                string sensspecFilename = fileroot+ tag + ".sensspec";
                ofstream sensFile;
                util.openOutputFile(sensspecFilename,    sensFile);
                outputNames.push_back(sensspecFilename); outputTypes["sensspec"].push_back(sensspecFilename);

                //evaluate results
                bestListFileName = compareSensSpec(matrix, metric, sensFile);

                delete matrix;

            }else { //reference with accnos file or reference list file assigning references

                unordered_set<string> refNames; vector<string> refLabels; vector< vector<string> > otus;

                if (accnosfile != "") { //use accnos file to assign references

                    m->mothurOut("\nUsing sequences from " + accnosfile + " as reference sequences\n");

                    refNames = util.readAccnos(accnosfile);

                }else if (createAccnos) { //assign references based on reflist parameter

                    m->mothurOut("\nUsing OTUs from " + reflistfile + " as reference OTUs\n");

                    InputData input(reflistfile, "list", nullVector);
                    set<string> processedLabels, userLabels;
                    string lastLabel = "";

                    ListVector* reflist = util.getNextList(input, true, userLabels, processedLabels, lastLabel);

                    refLabels = reflist->getLabels();
                    for (int i = 0; i < refLabels.size(); i++) { refLabels[i] = "Ref_" + refLabels[i];  }

                    refNames = util.getSetFromList(reflist, otus); delete reflist;
                }

                //distfile, distFormat, dupsFile, dupsFormat, cutoff, accnos containing refseq name
                OptiData* matrix = new OptiRefMatrix(distfile, "column", dupsFile, nameOrCount, cutoff, refNames);

                //fit seqs
                ListVector* list = runUserRefOptiCluster(matrix, metric, counts, outputName, refLabels, otus);

                ofstream listFile; string listFileName = fileroot+ tag + ".list";
                util.openOutputFile(listFileName,    listFile);

                if(countfile != "") { list->print(listFile, counts); }
                else { list->print(listFile); }
                listFile.close();

                listFiles.push_back(listFileName);
                bestListFileName = listFileName;

                delete list; delete matrix;
            }
        }else { //reference with files containing reference seqs

            createReferenceNameCount(); //creates reference name or count file if needed

            string distanceFile = calcDists();  //calc distance matrix for fasta file and distances between fasta file and reffasta file

            if (outputdir == "") { outputdir += util.hasPath(distanceFile); }
            fileroot = outputdir + util.getRootName(util.getSimpleName(distanceFile));

            map<string, string> variables;
            variables["[filename]"] = fileroot;
            variables["[clustertag]"] = "optifit_" + metric->getName();
            outputName = getOutputFileName("steps", variables);

            m->mothurOut("\nUsing OTUs from " + reflistfile + " as reference OTUs\n");

            //calc sens.spec values for reference
            InputData input(reflistfile, "list", nullVector);
            ListVector* list = input.getListVector();

            //add tag to OTULabels to indicate the reference
            vector<string> refListLabels = list->getLabels();
            for (int i = 0; i < refListLabels.size(); i++) { refListLabels[i] = "Ref_" + refListLabels[i];  }
            list->setLabels(refListLabels);

            string refDupsFile = refcountfile;
            if (refNameOrCount == "name") { refDupsFile = refnamefile; }

            OptiData* matrix = new OptiRefMatrix(refdistfile, refDupsFile, refNameOrCount, refformat, cutoff, distfile, dupsFile, nameOrCount, "column", comboDistFile, "column");

            listFile = runRefOptiCluster(matrix, metric, list, counts, outputName);
            listFiles.push_back(listFile);

            bestListFileName = listFile;

            delete matrix;
        }
        delete metric;

        if (m->getControl_pressed()) {     for (int j = 0; j < outputNames.size(); j++) { util.mothurRemove(outputNames[j]); }  return 0; }

        outputNames.push_back(outputName); outputTypes["steps"].push_back(outputName);
        outputNames.push_back(bestListFileName); outputTypes["list"].push_back(bestListFileName);

        if (m->getControl_pressed()) {     for (int j = 0; j < outputNames.size(); j++) { util.mothurRemove(outputNames[j]); }  return 0; }

        m->mothurOut("It took " + toString(time(nullptr) - estart) + " seconds to fit sequences to reference OTUs.\n");

        //set list file as new current listfile
        string currentName = "";
        itTypes = outputTypes.find("list");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setListFile(currentName); }
        }

        itTypes = outputTypes.find("accnos");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setAccnosFile(currentName); }
        }

        m->mothurOut("\nOutput File Names: \n");
        for (int i = 0; i < outputNames.size(); i++) {    m->mothurOut(outputNames[i]+"\n");     }
        m->mothurOutEndLine();

        return 0;
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

        int binNumber = 0;
        int placeHolderIndex = -1;
        for (long long i = 0; i < translatedBins.size(); i++) {
            binLabels[binNumber] = bls[i];
            bins.push_back(translatedBins[i]);
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
        int numRefBins = translatedBins.size();
        numFitSingletons = 0;
        //put every fit seq in own bin
        for (long long i = 0; i < randomizeSeqs.size(); i++) {
            std::vector<long long> thisBin;
            thisBin.push_back(randomizeSeqs[i]);
            bins.push_back(thisBin);
            seqBin[randomizeSeqs[i]] = numRefBins+i;

            long long numCloseSeqs = (matrix->getNumFitClose(randomizeSeqs[i])); //does not include self
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

        double comboValue = metric->getValue(combotruePositives, combotrueNegatives, combofalsePositives, combofalseNegatives);

        //add insert location
        seqBin[bins.size()] = -1;
        insertLocation = bins.size();
        std::vector<long long> temp;
        bins.push_back(temp);

        if (randomize) { Utils::mothurRandomShuffle(randomizeSeqs); }

        value = comboValue;

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
        ListVector* list = new ListVector();
        const auto* singleton = new ListVector(matrix->getListSingle());

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

        Rcpp::Rcout << ("\nFitted " + std::to_string(numListSeqs) + " sequences to " + std::to_string(newBins.size()) + " existing OTUs.\n");

        if (!closed) { //cluster the unfitted seqs separately
            Rcpp::Rcout << (std::to_string(numUnFitted) + " sequences were unable to be fitted existing OTUs, excluding singletons.\n");

            Rcpp::Rcout << ("\n**************** Clustering the unfitted sequences ****************\n");

            OptiRefMatrix* unFittedMatrix = matrix->extractMatrixSubset(unFitted);
            //matrix->extractMatrixSubset(unFitted)

            ListVector* unfittedList = clusterUnfitted(unFittedMatrix, label); //unfittedList includes unfitted singletons

            if (unfittedList != nullptr) {

                Rcpp::Rcout << ("The unfitted sequences clustered into " + std::to_string(unfittedList->getNumBins()) + " new OTUs.\n"); //+unFittedMatrix->getNumSingletons()+ matrix->getNumFitSingletons()

                for (int i = 0; i < unfittedList->getNumBins(); i++) {
                    std::string bin = unfittedList->get(i);
                    if (bin != "") { list->push_back(unfittedList->get(i)); }
                }
                delete unfittedList;
            }
            delete unFittedMatrix;

            Rcpp::Rcout << ("\n*******************************************************************\n\n");

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
            Rcpp::Rcout << ("\nSequences that were unable to be fitted existing OTUs will be listed in the *.optifit_scrap.accnos file.\n");
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
            Rcpp::Rcout << ("\nFitted all " + std::to_string(list->getNumSeqs()) + " sequences to existing OTUs. \n");
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
ListVector* OptiFitCluster::clusterUnfitted(OptiRefMatrix* unfittedMatrix, std::string label) {

    ListVector* list = nullptr;
    auto* opti = new OptiMatrix();
    opti->MoveData(unfittedMatrix);
    OptiCluster cluster(opti, metric, opti->GetCutoff(), 0);
    int iters = 0;
    double listVectorMetric = 0; //worst state
    double delta = 1;

    cluster.initialize(listVectorMetric, true, "singleton");

    long long numBins = cluster.getNumBins();
    Rcpp::Rcout << ("\n\niter\ttime\tlabel\tnum_otus\tcutoff\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\n");

    double tp, tn, fp, fn;
    std::vector<double> results = cluster.getStats(tp, tn, fp, fn);
    Rcpp::Rcout << ("0\t0\t" + label + "\t" + std::to_string(numBins) + "\t"+ label + "\t" + std::to_string(tp) + "\t" + std::to_string(tn) + "\t" + std::to_string(fp) + "\t" + std::to_string(fn) + "\t");

    for (const double result : results) { Rcpp::Rcout << (std::to_string(result) + "\t");  }

    while ((delta > 0.0001) && (iters < 100)) {

        const long start = time(nullptr);

        const double oldMetric = listVectorMetric;

        cluster.update(listVectorMetric);

        delta = abs(oldMetric - listVectorMetric);
        iters++;

        results = cluster.getStats(tp, tn, fp, fn);
        numBins = cluster.getNumBins();

        Rcpp::Rcout << (std::to_string(iters) + "\t" + std::to_string(time(nullptr) - start) + "\t" + label + "\t" + std::to_string(numBins) + "\t" + label + "\t"+ std::to_string(tp) + "\t" + std::to_string(tn) + "\t" + std::to_string(fp) + "\t" + std::to_string(fn) + "\t");

        for (const double result : results) { Rcpp::Rcout << (std::to_string(result) + "\t");  }

    }
    list = new ListVector(cluster.getList());
    list->setLabel(label);

    return list;
}

/***********************************************************************/
long long OptiFitCluster::getNumBins() {
        long long singletn = 0;

        singletn = matrix->getNumSingletons();

        for (int i = 0; i < bins.size(); i++) { if (bins[i].size() != 0) { singletn++; } }

        return singletn;
}
/***********************************************************************/
long long OptiFitCluster::getNumFitBins() {
        ListVector* list = getFittedList("", false);

        int numBins = 0;
        if (list != nullptr) {
            numBins = list->getNumBins();
            delete list;
        }

        return numBins;
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