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
    ListVectorPair(ListVector * const list_vector, std::string label)
        : listVector(list_vector),
          label(std::move(label)) {
    }
    ListVector* listVector;
    std::string label;
};
class ClusterExport {
public:

    ClusterExport() = default;
    virtual void SetListVector(ListVector& vector, const std::string& label)  {
        listVector = &vector;
        largestLabel = label;
    };
    virtual ListVectorPair GetListVector() const {return {listVector, largestLabel};}
    virtual int GetLargestBinSize() const {
        return listVector->getNumBins();
    }
    virtual ~ClusterExport() {
        delete listVector;
    }
protected:
    ListVector* listVector = nullptr;
    std::string largestLabel;
};

#endif //CLUSTEREXPORT_H
