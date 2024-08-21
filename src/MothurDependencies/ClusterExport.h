//
// Created by gregj on 8/20/2024.
//

#ifndef CLUSTEREXPORT_H
#define CLUSTEREXPORT_H
#include <string>

#include "RAbundVector.h"

class ClusterExport {
public:
    virtual ~ClusterExport() = default;
    ClusterExport() = default;
    virtual std::string Print() const = 0;
    virtual std::string GetLabel(int) const = 0;
    virtual size_t GetAbundanceCount() = 0;
};

#endif //CLUSTEREXPORT_H
