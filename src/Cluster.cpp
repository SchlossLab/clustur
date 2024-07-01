//
// Created by Gregory Johnson on 6/14/24.
//

#include "MothurDependencies/Cluster.h"
/*
 *  cluster.cpp
 *
 *
 *  Created by Pat Schloss on 8/14/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */


/***********************************************************************/

Cluster::Cluster(RAbundVector* rav, ListVector* lv, SparseDistanceMatrix* dm, float c, std::string f, float cs) :
rabund(rav), list(lv), dMatrix(dm), method(f), adjust(cs)
{
        mapWanted = false;  //set to true by mgcluster to speed up overlap merge

        //save so you can modify as it changes in average neighbor
        cutoff = c;
}
/***********************************************************************/
void Cluster::clusterBins(){

 		rabund->set(smallCol, rabund->get(smallRow)+rabund->get(smallCol));
		rabund->set(smallRow, 0);
		rabund->setLabel(std::to_string(smallDist));

}
/***********************************************************************/

void Cluster::clusterNames(){

		if (mapWanted) {  updateMap();  }

		list->set(smallCol, list->get(smallRow)+','+list->get(smallCol));
		list->set(smallRow, "");
		list->setLabel(std::to_string(smallDist));
    // std::ofstream stream;
    // std::cout << list->print(stream);

}
/***********************************************************************/
bool Cluster::update(double& cutOFF){
        smallCol = dMatrix->getSmallestCell(smallRow);
        if(smallCol == -1)
            return false;
        nColCells = dMatrix->seqVec[smallCol].size();
        nRowCells = dMatrix->seqVec[smallRow].size();
        // std::cout << "rowCells: " << nRowCells << std::endl;
		std::vector<int> foundCol(nColCells, 0);

		int search;
		bool changed = false;

		for (int i=nRowCells-1;i>=0;i--) {  //matrix indexes sorted from largest to smallest, so start at smallest index


			//if you are not the smallCell
			if (dMatrix->seqVec[smallRow][i].index != smallCol) {
                search = dMatrix->seqVec[smallRow][i].index;

				bool merged = false;
				for (int j=0;j<nColCells;j++) {  //go through each distance the smallCol has looking for matching distance to find

					if (dMatrix->seqVec[smallCol][j].index != smallRow) {  //if you are not the smallest distance
						if (dMatrix->seqVec[smallCol][j].index == search) {  //we found a distance for the merge
							foundCol[j] = 1;
							merged = true;
							changed = updateDistance(dMatrix->seqVec[smallCol][j], dMatrix->seqVec[smallRow][i]);
                            dMatrix->updateCellCompliment(smallCol, j);
							break;
						}else if (dMatrix->seqVec[smallCol][j].index < search) { //we don't have a distance for this cell
                            if (!util.isEqual(adjust, -1)) { //adjust
                                merged = true;
                                PDistCell value(search, adjust); //create a distance for the missing value
                                int location = dMatrix->addCellSorted(smallCol, value);
                                changed = updateDistance(dMatrix->seqVec[smallCol][location], dMatrix->seqVec[smallRow][i]);
                                dMatrix->updateCellCompliment(smallCol, location);
                                nColCells++;
                                foundCol.push_back(0); //add a new found column
                                //adjust value
                                for (int k = foundCol.size()-1; k > location; k--) { foundCol[k] = foundCol[k-1]; }
                                foundCol[location] = 1;
                            }
                            j+=nColCells;  //jump out of loop and remove cell below
                        }
                    }
				}
				//if not merged it you need it for warning
				if ((!merged) && (method == "average" || method == "weighted")) {   if (cutOFF > dMatrix->seqVec[smallRow][i].dist) {   cutOFF = dMatrix->seqVec[smallRow][i].dist; } }
                if ((method == "nearest") && (!merged)) { //you are a row dist without a column dist, add you as a column dist
                    PDistCell value(search, dMatrix->seqVec[smallRow][i].dist); //create a distance for the missing value
                    int location = dMatrix->addCellSorted(smallCol, value); nColCells++;
                    foundCol.push_back(0); //add a new found column
                    //adjust value
                    for (int k = foundCol.size()-1; k > location; k--) { foundCol[k] = foundCol[k-1]; }
                    foundCol[location] = 1;
                }
                dMatrix->rmCell(smallRow, i);
			}
		}
		clusterBins();
		clusterNames();

        if (method == "nearest") {
            for (int i=nColCells-1;i>=0;i--) { //remove any unfound dists from merged column, need special case for nn, since unfound dists mean above the cutoff -> keep smaller dist in col
                if (foundCol[i] == 0) {  //not found
                    if (dMatrix->seqVec[smallCol][i].index == smallRow) { //you are smallest distance
                        dMatrix->rmCell(smallCol, i);
                        break;
                    }
                }
            }
        }else {
            for (int i=nColCells-1;i>=0;i--) { //remove any unfound dists from merged column, need special case for nn, since unfound dists mean above the cutoff -> keep smaller dist in col
                if (foundCol[i] == 0) {  //not found
                    if (!util.isEqual(adjust, -1)) { //adjust
                        PDistCell value(smallCol, adjust); //create a distance for the missing value
                        changed = updateDistance(dMatrix->seqVec[smallCol][i], value);
                        dMatrix->updateCellCompliment(smallCol, i);
                    }else {
                        if (method == "average" || method == "weighted") {
                            if (dMatrix->seqVec[smallCol][i].index != smallRow) { //if you are not hte smallest distance
                                if (cutOFF > dMatrix->seqVec[smallCol][i].dist) {   cutOFF = dMatrix->seqVec[smallCol][i].dist;   }
                            }
                        }
                    }
                    dMatrix->rmCell(smallCol, i);
                }
            }
        }
        //dMatrix->print();
        return changed;
}
/***********************************************************************/
void Cluster::setMapWanted(bool f)  {
		mapWanted = f;

        //initialize map
		for (int k = 0; k < list->getNumBins(); k++) {

            std::string names = list->get(k);

            //parse bin
            std::string individual = "";
            int binNameslength = names.size();
            for(int j=0;j<binNameslength;j++){
                if(names[j] == ','){
                    seq2Bin[individual] = k;
                    individual = "";
                }
                else{  individual += names[j];  }
            }
            //get last name
            seq2Bin[individual] = k;
		}

}
/***********************************************************************/
void Cluster::updateMap() {
		//update location of seqs in smallRow since they move to smallCol now
		std::string names = list->get(smallRow);

        std::string individual = "";
        int binNameslength = names.size();
        for(int j=0;j<binNameslength;j++){
            if(names[j] == ','){
                seq2Bin[individual] = smallCol;
                individual = "";
            }
            else{  individual += names[j];  }
        }
        //get last name
        seq2Bin[individual] = smallCol;


}
/***********************************************************************/



