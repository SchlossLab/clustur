//
// Created by Gregory Johnson on 6/25/26.
//

#ifndef REFACTOR_CLUSTERTYPE_H
#define REFACTOR_CLUSTERTYPE_H
#include <string>
#include "../DataStructures/PDistCell.h"
class ClusterType {
public:
    ClusterType() = default;
    virtual ~ClusterType() = default;
    virtual bool updateDistance(PDistCell& colCell, PDistCell& rowCell) = 0;
    virtual std::string getTag() = 0;
};
#endif //REFACTOR_CLUSTERTYPE_H