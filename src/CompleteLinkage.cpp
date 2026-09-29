//
// Created by Gregory Johnson on 6/17/24.
//

#include "Clusters/CompleteLinkage.h"

/***********************************************************************/

CompleteLinkage::CompleteLinkage(RAbundVector* rav, ListVector* lv, SparseDistanceMatrix* dm, const float c,
    const std::string &s, const float a) :
    Cluster(rav, lv, dm, c, s, a)
{}

/***********************************************************************/
//This function returns the tag of the method.
std::string CompleteLinkage::getTag() {
    return("fn");
}

ClusterExport* CompleteLinkage::Execute() {
    return ExecuteCluster();
}

void CompleteLinkage::SetClusterParameters(const ClusterParameters &parameters) {
    const std::unordered_map<std::string, std::string>& params = parameters.GetClusterParameters();
    precision = std::stoi(params.at("precision"));
}


/***********************************************************************/
//This function updates the distance based on the furthest neighbor method.
bool CompleteLinkage::updateDistance(PDistCell& colCell, PDistCell& rowCell) {

    bool changed = false;
    if (colCell.dist < rowCell.dist) {
        colCell.dist = rowCell.dist;
        changed = true;
    }
    return changed;

}

/***********************************************************************/
