//
// Created by Gregory Johnson on 9/2/26.
//

#ifndef REFACTOR_CLUSTERSPLIT_H
#define REFACTOR_CLUSTERSPLIT_H
#include "ClusterMethod.h"


class ClusterSplit : ClusterMethod {
public:
    ClusterSplit() = default;
    ~ClusterSplit() override;
    ClusterExport * Execute() override;
};


#endif //REFACTOR_CLUSTERSPLIT_H