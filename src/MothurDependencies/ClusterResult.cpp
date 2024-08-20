//
// Created by Gregory Johnson on 8/19/24.
//

#include "ClusterResult.h"


std::string ClusterResult::ResultsToString(const ClusterExport& exportedData) const {
    return exportedData.Print();
}
