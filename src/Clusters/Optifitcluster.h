
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

class OptiFitCluster : public ClusterMethod {

public:

    OptiFitCluster(OptiRefMatrix* mt, ClusterMetric* met, long long ns);
    ~OptiFitCluster() override = default;

    int initialize(double& value, bool randomize, std::vector<std::vector<std::string>>& existingBins,
        const std::vector<std::string>&, const std::string&, bool);
    bool update(double&); //returns whether list changed and MCC
    [[nodiscard]] std::string getTag() const { std::string tag = "optifit_" + metric->getName(); return tag; }
    long long getNumBins();
    long long getNumFitBins();

    std::vector<double>  getStats( double&, double&,  double&,  double&);  //combo stats
    std::vector<double>  getFitStats( double&,  double&,  double&,  double&); //fitted seqs stats

    ListVector* getList();
    ListVector* getFittedList(std::string, bool);
    std::set<std::string> getUnfittedNames() { return unfittedNames; }

protected:
    ClusterMetric* metric;
    OptiRefMatrix* matrix;

    std::map<long long, long long> seqBin; //sequence# -> bin#
    std::vector<long long> randomizeSeqs;
    std::vector< std::vector<long long> > bins; //bin[0] -> seqs in bin[0]
    std::map<long long, std::string> binLabels; //for fitting - maps binNumber to existing reference label
    long long maxRefBinNumber;
    bool closed, denovo;
    std::set<std::string> unfittedNames;

    double fittruePositives, fittrueNegatives, fitfalsePositives, fitfalseNegatives, combotruePositives, combotrueNegatives, combofalsePositives, combofalseNegatives;
    long long  numFitSeqs, insertLocation, numFitSingletons;
    long long  numComboSeqs, numComboSingletons;

    int findInsert();
    std::vector<double> getCloseFarCounts(long long seq, long long newBin);
    std::vector<double> getCloseFarFitCounts(long long seq, long long newBin);
    ListVector* clusterUnfitted(OptiRefMatrix*, std::string);

};


#endif //REFACTOR_OPTIFITCLUSTER_H