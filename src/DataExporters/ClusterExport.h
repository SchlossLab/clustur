//
// Created by Gregory Johnson on 6/9/26.
//

#ifndef REFACTOR_CLUSTEREXPORT_H
#define REFACTOR_CLUSTEREXPORT_H


#include <string>
#include <utility>

#include "../DataStructures/ListVector.h"
#include "../DataStructures/RAbundVector.h"
struct ListVectorPair {
    ListVector listVector;
    std::string label;
};
class ClusterExport {
public:

    ClusterExport() = default;
     void SetListVector(const ListVector& vector, const std::string& label)  {
        listVector = vector;
        largestLabel = label;
    };
    [[nodiscard]] ListVectorPair GetListVector() const { return ListVectorPair{listVector, largestLabel};}
    [[nodiscard]] int GetLargestBinSize() const {
        return listVector.getNumBins();
    }
    virtual ~ClusterExport() = default;
protected:
    ListVector listVector{};
    std::string largestLabel;
};


#endif //REFACTOR_CLUSTEREXPORT_H