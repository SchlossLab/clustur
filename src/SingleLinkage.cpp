//
// Created by Gregory Johnson on 6/17/24.
//

#include "Clusters/SingleLinkage.h"

/***********************************************************************/

SingleLinkage::SingleLinkage(RAbundVector* rav, ListVector* lv, SparseDistanceMatrix* dm, const float c,
    const std::string& s, const float a) :
Cluster(rav, lv, dm, c, s, a)
{}


/***********************************************************************/
//This function returns the tag of the method.
std::string SingleLinkage::getTag() {
    return("nn");
}

ClusterExport* SingleLinkage::Execute() {
    return ExecuteCluster();
}

void SingleLinkage::SetClusterParameters(const ClusterParameters &parameters) {
    const std::unordered_map<std::string, std::string>& params = parameters.GetClusterParameters();
    precision = std::stoi(params.at("precision"));
}

/***********************************************************************/
//This function updates the distance based on the nearest neighbor method.
bool SingleLinkage::updateDistance(PDistCell& colCell, PDistCell& rowCell) {
    constexpr bool changed = false;
    if (colCell.dist > rowCell.dist) {
        colCell.dist = rowCell.dist;
    }
    return changed;

}
/***********************************************************************/
