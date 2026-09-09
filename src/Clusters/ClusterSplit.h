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

class ClusterSplit : ClusterMethod {
public:
    ClusterSplit(const FastaDatabase& fastaDatabase, const std::vector<TaxonomyData>& taxaData,
        PairwiseDistanceCalculator* calculator, double cutoff, int taxonomyCutoff);
    ~ClusterSplit() override;
    ClusterExport * Execute() override;

    void cluster(ClusterData *params);

private:
    FastaDatabase fastaData;
    PairwiseDistanceCalculator* calculator;
    std::vector<TaxonomyData> taxaData;
    int taxonomyCutoff;
    double cutoff;

};


#endif //REFACTOR_CLUSTERSPLIT_H