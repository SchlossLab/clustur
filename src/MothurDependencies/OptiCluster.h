//
// Created by Gregory Johnson on 3/29/24.
//

#ifndef OPTICLUSTER_H
#define OPTICLUSTER_H
#include "OptiData.h"
#include "Utils.h"
#include "ClusterMetric.h"
#include "ListVector.h"
#include <map>
#include <unordered_map>
#include <string>
#include <set>
#include <vector>

class OptiCluster {
public:

    OptiCluster(OptiData* mt, ClusterMetric* met, long long ns);
    ~OptiCluster() = default;

    std::string getTag() { std::string tag = "opti_" + metric->getName(); return tag; }
    long long getNumBins() const;
    int initialize(double&, bool, const std::string&);  //randomize and place in "best" OTUs

    bool update(double&); //returns whether list changed and MCC
    std::vector<double> getStats( double&,  double&,  double&,  double&) const;
    std::vector<double> getCloseFarFitCounts(long long seq, long long newBin) const;
    std::vector<double> getCloseFarCounts(long long seq, long long newBin) const;
    ListVector* getList() const;

protected:
    OptiData* matrix;
    ClusterMetric* metric;
    std::vector<int> randomizeSeqs;
    std::vector< std::vector<long long> > bins; //bin[0] -> seqs in bin[0]
    std::map<long long, std::string> binLabels; //for fitting - maps binNumber to existing reference label
    std::unordered_map<long long, long long> seqBin; //sequence# -> bin#

    long long numSeqs, insertLocation, numSingletons;
    double fittruePositives, fittrueNegatives, fitfalsePositives, fitfalseNegatives, combotruePositives, combotrueNegatives, combofalsePositives, combofalseNegatives;
    long long  numFitSeqs, numFitSingletons;
    long long  numComboSeqs, numComboSingletons;
    double truePositives, trueNegatives, falsePositives, falseNegatives;
    long long findInsert() const;
};



#endif //OPTICLUSTER_H
