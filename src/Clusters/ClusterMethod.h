//
// Created by Gregory Johnson on 6/9/26.
//

#ifndef REFACTOR_CLUSTERMETHOD_H
#define REFACTOR_CLUSTERMETHOD_H
#include "../DataExporters/ClusterExport.h"
#include "../Adapters/DataFrameAdapter.h"

class ClusterMethod {
public:
    ClusterMethod() = default;
    virtual ~ClusterMethod() = default;
    virtual ClusterExport* Execute() = 0;
    Rcpp::DataFrame GetSensitivityData() const {return DataFrameAdapter::UnorderedMapToDataFrame(dataframeMapSensMetrics);}
    Rcpp::DataFrame GetClusterMetrics() const {return DataFrameAdapter::UnorderedMapToDataFrame(dataframeMapClusterMetrics);}
protected:
    std::unordered_map<std::string, std::vector<std::string>> dataframeMapSensMetrics;
    std::unordered_map<std::string, std::vector<std::string>> dataframeMapClusterMetrics;
};
#endif //REFACTOR_CLUSTERMETHOD_H