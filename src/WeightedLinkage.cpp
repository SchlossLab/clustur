//
// Created by Gregory Johnson on 6/17/24.
//

#include <utility>

#include "MothurDependencies/WeightedLinkage.h"
/* This class implements the WPGMA, weighted average neighbor clustering algorithm */

/***********************************************************************/

WeightedLinkage::WeightedLinkage(RAbundVector* rav, ListVector* lv, SparseDistanceMatrix* dm, const float c, std::string s, const float a) :
    Cluster(rav, lv, dm, c, std::move(s), a)
{
    saveRow = -1;
    saveCol = -1;
}


/***********************************************************************/
//This function returns the tag of the method.
std::string WeightedLinkage::getTag() {
    return("wn");
}


/***********************************************************************/
//This function updates the distance based on the average linkage method.
bool WeightedLinkage::updateDistance(PDistCell& colCell, PDistCell& rowCell) {

    if ((saveRow != smallRow) || (saveCol != smallCol)) {
        saveRow = smallRow;
        saveCol = smallCol;
    }

    colCell.dist = (colCell.dist + rowCell.dist) / 2;

    return(true);

}

/***********************************************************************/


/***********************************************************************/
