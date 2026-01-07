//
// Created by gregj on 8/20/2024.
//

#ifndef CLUSTEREXPORT_H
#define CLUSTEREXPORT_H
#include <string>
#include <utility>

#include "ListVector.h"
#include "RAbundVector.h"
struct ListVectorPair {
    ListVector listVector;
    std::string label;
};
class ClusterExport {
public:

    ClusterExport() = default;
    virtual void SetListVector(const ListVector& vector, const std::string& label)  {
        listVector = vector;
        largestLabel = label;
    };
    virtual ListVectorPair GetListVector() const {return ListVectorPair{listVector, largestLabel};}
    virtual int GetLargestBinSize() const {
        return listVector.getNumBins();
    }
    virtual ~ClusterExport() = default;
protected:
    ListVector listVector;
    std::string largestLabel;
};

#endif //CLUSTEREXPORT_H
