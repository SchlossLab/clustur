//
// Created by Gregory Johnson on 6/17/24.
//

#include <utility>

#include "MothurDependencies/AverageLinkage.h"

/* This class implements the average UPGMA, average neighbor clustering algorithm */

/***********************************************************************/

AverageLinkage::AverageLinkage(RAbundVector* rav, ListVector* lv, SparseDistanceMatrix* dm,
const float c, std::string s, const float a) :
Cluster(rav, lv, dm, c, std::move(s), a),saveRow(-1), saveCol(-1), rowBin(-1), colBin(-1),totalBin(-1)
{}


/***********************************************************************/
//This function returns the tag of the method.
std::string AverageLinkage::getTag() {
    return("an");
}


/***********************************************************************/
//This function updates the distance based on the average linkage method.
bool AverageLinkage::updateDistance(PDistCell& colCell, PDistCell& rowCell) {
    if ((saveRow != smallRow) || (saveCol != smallCol)) {
        rowBin = rabund->get(static_cast<int>(smallRow));
        colBin = rabund->get(static_cast<int>(smallCol));
        totalBin = rowBin + colBin;
        saveRow = smallRow;
        saveCol = smallCol;
    }

    colCell.dist = (static_cast<float>(colBin) * colCell.dist + static_cast<float>(rowBin) * rowCell.dist) /
        static_cast<float>(totalBin);

    return(true);
}