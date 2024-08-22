//
// Created by gregj on 8/20/2024.
//

#ifndef CLUSTEREXPORT_H
#define CLUSTEREXPORT_H
#include <string>

#include "ListVector.h"
#include "RAbundVector.h"

class ClusterExport {
public:
    virtual ~ClusterExport() = default;
    ClusterExport() = default;
    virtual std::string Print() const = 0;
    virtual std::string GetLabel(int) const = 0;
    virtual void AddToList(ListVector& vector)  {
        listVectors.emplace_back(&vector);
    };
    virtual std::vector<ListVector*> GetListVectors() const {return listVectors;}
protected:
    std::vector<ListVector*> listVectors;

};

#endif //CLUSTEREXPORT_H
