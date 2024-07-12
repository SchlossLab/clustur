//
// Created by Gregory Johnson on 6/17/24.
//

#include "MothurDependencies/CompleteLinkage.h"

/***********************************************************************/

CompleteLinkage::CompleteLinkage(RAbundVector* rav, ListVector* lv, SparseDistanceMatrix* dm, float c, std::string s, float a) :
    Cluster(rav, lv, dm, c, s, a)
{}

/***********************************************************************/
//This function returns the tag of the method.
std::string CompleteLinkage::getTag() {
    return("fn");
}


/***********************************************************************/
//This function updates the distance based on the furthest neighbor method.
bool CompleteLinkage::updateDistance(PDistCell& colCell, PDistCell& rowCell) {

    bool changed = false;
    if (colCell.dist < rowCell.dist) {
        colCell.dist = rowCell.dist;
        changed = true;
    }
    return(changed);

}

/***********************************************************************/
