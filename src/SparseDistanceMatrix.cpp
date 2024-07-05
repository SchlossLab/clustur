//
// Created by Gregory Johnson on 6/17/24.
//

#include "MothurDependencies/SparseDistanceMatrix.h"
#include <iostream>
#include "MothurDependencies/sparsedistancematrix.h"


/***********************************************************************/

SparseDistanceMatrix::SparseDistanceMatrix() : numNodes(0), smallDist(MOTHURMAX){ sorted=false; aboveCutoff = MOTHURMAX;}

/***********************************************************************/

int SparseDistanceMatrix::getNNodes(){
	return numNodes; 
}
/***********************************************************************/

void SparseDistanceMatrix::clear(){
    for (int i = 0; i < seqVec.size(); i++) {  seqVec[i].clear();  }
    seqVec.clear();
}

/***********************************************************************/

float SparseDistanceMatrix::getSmallDist(){
	return smallDist;
}

bool SparseDistanceMatrix::heapComparator(const PDistCell &a, const PDistCell &b) {
    return a.dist < b.dist;
}

/***********************************************************************/

int SparseDistanceMatrix::updateCellCompliment(unsigned long row, unsigned long col){

        
    unsigned long vrow = seqVec[row][col].index;
    unsigned long vcol = 0;

    //find the columns entry for this cell as well
    for (int i = 0; i < seqVec[vrow].size(); i++) {
        if (seqVec[vrow][i].index == row) { vcol = i;  break; }
    }

    seqVec[vrow][vcol].dist = seqVec[row][col].dist;

    return 0;

}
/***********************************************************************/

int SparseDistanceMatrix::rmCell(unsigned long row, unsigned long col){
    numNodes-=2;

    unsigned long vrow = seqVec[row][col].index;
    unsigned long vcol = 0;

    //find the columns entry for this cell as well
    for (int i = 0; i < seqVec[vrow].size(); i++) {  if (seqVec[vrow][i].index == row) { vcol = i;  break; }  }
    seqVec[vrow].erase(seqVec[vrow].begin()+vcol);
    seqVec[row].erase(seqVec[row].begin()+col);


	return(0);
}
/***********************************************************************/
void SparseDistanceMatrix::addCell(const unsigned long row, const PDistCell cell){

    numNodes+=2;

    if(cell.dist < smallDist){ smallDist = cell.dist;}
    seqVec[row].push_back(cell);
    if(row == cell.index) return;
    const PDistCell temp(row, cell.dist);
    seqVec[cell.index].push_back(temp);

    // It pushes itself back in the row, and in the column to create the sparse matrix
}
/***********************************************************************/
int SparseDistanceMatrix::addCellSorted(unsigned long row, PDistCell cell){

	numNodes+=2;
	if(cell.dist < smallDist){ smallDist = cell.dist; }

    seqVec[row].push_back(cell);
    PDistCell temp(row, cell.dist);
    seqVec[cell.index].push_back(temp);

    sortSeqVec(row);
    sortSeqVec(cell.index);

    int location = -1; //find location of new cell when sorted
    for (int i = 0; i < seqVec[row].size(); i++) {  if (seqVec[row][i].index == cell.index) { location = i; break; } }

    return location;
}

/***********************************************************************/

unsigned long SparseDistanceMatrix::getSmallestCell(unsigned long& row){
    if (!sorted) { sortSeqVec(); sorted = true; }

    //print();

    std::vector<PDistCellMin> mins;
    smallDist = MOTHURMAX;
    for (int i = 0; i < seqVec.size(); i++) {
        for (int j = 0; j < seqVec[i].size(); j++) {

            if (i < seqVec[i][j].index) {
                const float dist = seqVec[i][j].dist;
                if(dist < smallDist){  //found a new smallest distance
                    mins.clear();
                    smallDist = dist;
                    PDistCellMin temp(i, seqVec[i][j].index);
                    mins.emplace_back(temp);
                }
                else if(util.isEqual(dist, smallDist)){  //if a subsequent distance is the same as mins distance add the new iterator to the mins vector
                    PDistCellMin temp(i, seqVec[i][j].index);
                    mins.emplace_back(temp);
                }
            }else { j+=seqVec[i].size(); } //stop looking
		}
	}
    if(mins.empty())
        return -1;
    // TestHelper::Print(std::to_string(smallDist) + "\n");
    const unsigned long num = util.getRandomNumber(static_cast<int>(mins.size() - 1));
    row = mins[num].row;
    const unsigned long col = mins[num].col;
	return col;

}
/***********************************************************************/

void SparseDistanceMatrix::print(){

    std::cout << std::endl;
    //saves time in getSmallestCell, by making it so you dont search the repeats
    for (int i = 0; i < seqVec.size(); i++) {
        std::cout << i << '\t';
        for (int j = 0; j < seqVec[i].size(); j++) {   std::cout <<  seqVec[i][j].index << '\t' ;  }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
/***********************************************************************/

int SparseDistanceMatrix::sortSeqVec(){
        
        //saves time in getSmallestCell, by making it so you dont search the repeats
    for (int i = 0; i < seqVec.size(); i++) {  sort(seqVec[i].begin(), seqVec[i].end(), PDistCell::CompareIndexes); }

    return 0;
}
/***********************************************************************/

int SparseDistanceMatrix::sortSeqVec(int index){

    //saves time in getSmallestCell, by making it so you dont search the repeats
    sort(seqVec[index].begin(), seqVec[index].end(), PDistCell::CompareIndexes);

    return 0;

}
/***********************************************************************/

