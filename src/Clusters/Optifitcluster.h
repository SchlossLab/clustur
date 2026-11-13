
//
//  optifitcluster.hpp
//  Mothur
//
//  Created by Sarah Westcott on 5/10/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//
#ifndef REFACTOR_OPTIFITCLUSTER_H
#define REFACTOR_OPTIFITCLUSTER_H

#include "ClusterMethod.h"
#include "../Clusters/ClusterMetric.h"
#include "../DataStructures/OptiRefMatrix.h"
#include "../DataStructures/PDistCell.h"

/***********************************************************************/

class OptiFitCluster final : public ClusterMethod {

public:
    OptiFitCluster(OptiRefMatrix* mt, ClusterMetric* met,
        double cutoff, long long ns, int seed);
    OptiFitCluster(OptiRefMatrix* mt, ClusterMetric* met, const ListVector& refListVector,
        double cutoff, long long ns, int seed);

    void Reset();

    ~OptiFitCluster() override = default;
    ClusterExport* Execute() override;
    int initialize(double& value, bool randomize, std::vector<std::vector<std::string>>& existingBins,
        const std::vector<std::string>&, bool);
    bool update(double&); //returns whether list changed and MCC
    [[nodiscard]] std::string getTag() const { std::string tag = "optifit_" + metric->getName(); return tag; }
    long long getNumBins();
    long long getNumFitBins();

    std::vector<double>  getStats( double&, double&,  double&,  double&);  //combo stats
    std::vector<double>  getFitStats( double&,  double&,  double&,  double&) const; //fitted seqs stats

    [[nodiscard]] ListVector* getList() const;
    ListVector* getFittedList(const std::string&, bool);
    // std::set<std::string> getUnfittedNames() { return unfittedNames; }
    const std::vector<ListVector>& GetAllDenovoResults() {return denovoListResults;}
    const std::vector<std::unordered_map<std::string, std::vector<std::string>>>&
    GetAllSenSpecMetrics() {return senSpecDenovo;}
    std::vector<std::unordered_map<std::string, std::vector<std::string>>>
    GetAllClusterMetrics() {return clusterMetricsDenovo;}
    [[nodiscard]] size_t GetOptimalDenovoIndex() const {return optimalDenovoIndex;}

protected:
    ClusterMetric* metric;
    OptiRefMatrix* matrix;
    ListVector listVector;
    std::vector<std::pair<long long, long long>> seqBin; //sequence# -> bin#
    std::vector<long long> randomizeSeqs;
    std::vector< std::vector<long long> > bins; //bin[0] -> seqs in bin[0]
    std::map<long long, std::string> binLabels; //for fitting - maps binNumber to existing reference label
    long long maxRefBinNumber;
    bool closed, denovo;
    std::string method;
    std::vector<std::string> unfittedNames;
    double cutoff;
    double fittruePositives, fittrueNegatives, fitfalsePositives, fitfalseNegatives, combotruePositives, combotrueNegatives, combofalsePositives, combofalseNegatives;
    long long  numFitSeqs, insertLocation, numFitSingletons;
    long long  numComboSeqs, numComboSingletons;
    RandomNumberSitmo rng;
    int seed;

    [[nodiscard]] int findInsert() const;

    ClusterExport *runDenovoOptiCluster(std::map<std::string, int> &counts,
                                        const std::string& outStepFile);

    ClusterExport *runUserRefOptiCluster();

    ListVector clusterRefs(OptiData *&refsMatrix, ClusterMetric *&metric) const;

    ClusterExport *runRefOptiCluster(ListVector refList,
                                     std::map<std::string, int> &counts, const std::string& outStepFile);

    void runSensSpec(OptiRefMatrix *&matrix, ListVector *&list, std::string &sensSpecFile) const;

    void compareSensSpec(OptiRefMatrix *&matrix, ClusterMetric *&userMetric,
                                std::vector<ListVector> listVectors);


    [[nodiscard]] std::vector<double> getCloseFarCounts(long long seq, long long newBin) const;
    [[nodiscard]] std::vector<double> getCloseFarFitCounts(long long seq, long long newBin) const;
    ListVector clusterUnfitted(OptiData *unfittedMatrix) const;



public:
    void SetClusterParameters(const ClusterParameters &) override;

private:
    bool selfReference;
    bool printRef;
    int denovoIters;
    int maxIters;
    size_t optimalDenovoIndex;
    std::vector<ListVector> denovoListResults;
    std::vector<std::unordered_map<std::string, std::vector<std::string>>> senSpecDenovo;
    std::vector<std::unordered_map<std::string, std::vector<std::string>>> clusterMetricsDenovo;
};


#endif //REFACTOR_OPTIFITCLUSTER_H