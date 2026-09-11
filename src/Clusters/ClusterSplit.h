//
// Created by Gregory Johnson on 9/2/26.
//

#ifndef REFACTOR_CLUSTERSPLIT_H
#define REFACTOR_CLUSTERSPLIT_H
#include "ClusterMethod.h"
#include "../DataStructures/FastaDatabase.h"
#include "../DataStructures/OptiData.h"
#include "../DataStructures/TaxonomyData.h"
#include "../MothurDependencies/PairwiseDistanceCalculator.h"
#include "../DataStructures/ClusterData.h"
#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>

class ClusterSplit final : ClusterMethod {
public:
    ClusterSplit(FastaDatabase  fastaDatabase, const std::vector<TaxonomyData>& taxaData,
        PairwiseDistanceCalculator* calculator, ClusterParameters parameters, ClusterMetric* metric,
        double cutoff, int taxonomyCutoff);
    ~ClusterSplit() override = default;
    ClusterExport * Execute() override;

    ClusterExport *mergeLists(const std::vector<ClusterExport *> &exportedResults);

    std::vector<ClusterExport *> createProcesses(std::vector<OptiDataComponent> &distanceMatrices,
                                                 std::set<std::string> &labels, size_t processors);

    void cluster(ClusterData *params) const;

private:
    FastaDatabase fastaData;
    PairwiseDistanceCalculator* calculator;
    ClusterParameters clusterParameters;
    ClusterMetric* metric;
    std::vector<TaxonomyData> taxaData;
    int taxonomyCutoff;
    double cutoff;

};


#endif //REFACTOR_CLUSTERSPLIT_H