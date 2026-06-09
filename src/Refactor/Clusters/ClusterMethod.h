//
// Created by Gregory Johnson on 6/9/26.
//

#ifndef REFACTOR_CLUSTERMETHOD_H
#define REFACTOR_CLUSTERMETHOD_H
#include "DataExporters/ClusterExport.h"

class ClusterMethod {
public:
    ClusterMethod() = default;
    virtual ~ClusterMethod() = default;
    virtual ClusterExport Execute() = 0;
};
#endif //REFACTOR_CLUSTERMETHOD_H