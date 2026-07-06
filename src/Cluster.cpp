//
// Created by Gregory Johnson on 6/9/26.
//

#include "Clusters/Cluster.h"
#include "DataExporters/ClusterData.h"
#include <utility>

/*
 *  cluster.cpp
 *
 *
 *  Created by Pat Schloss on 8/14/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */


/***********************************************************************/

Cluster::Cluster(RAbundVector* rav, ListVector* lv, SparseDistanceMatrix* dm, const float c, std::string f,
                 const float cs) : rabund(rav), list(lv), dMatrix(dm), method(std::move(f)), adjust(cs) {
    mapWanted = false; //set to true by mgcluster to speed up overlap merge

    //save, so you can modify as it changes in average neighbor
    cutoff = c;
}

/***********************************************************************/
bool Cluster::clusterBins() {
    // if(list->size() <= 0)
    //     return false;
    const int sCol = static_cast<int>(smallCol);
    const int sRow = static_cast<int>(smallRow);
    rabund->set(sCol, rabund->get(sRow) + rabund->get(sCol));
    rabund->set(sRow, 0);
    rabund->setLabel(std::to_string(smallDist));
    return true;
}

/***********************************************************************/

bool Cluster::clusterNames() {
    // if(list->size() <= 0)
    //     return false;
    if (mapWanted) { updateMap(); }
    const int sCol = static_cast<int>(smallCol);
    const int sRow = static_cast<int>(smallRow);
    list->set(sCol, list->get(sRow) + ',' + list->get(sCol));
    list->set(sRow, "");
    list->setLabel(std::to_string(smallDist));
    return true;
}

ClusterExport* Cluster::ExecuteCluster() {
    double currentCutoff = cutoff;
    std::map<std::string, int> counts;
    float previousDist = 0.00000;
    float rndPreviousDist = 0.00000;
    ListVector oldList = *list;
    std::string clusterResult;
    double highestDistLabel =  -1;
    std::string binResults;
    std::ofstream listFile;
    auto* clusterData = new ClusterData("");
    // 0.158123 <= 0.1511
    while ((dMatrix->getSmallDist() <= currentCutoff) && (dMatrix->getNNodes() > 0)) {
        constexpr double precision = 100;
        //TODO We are getting values that are just barely grater than 0, we need to figure out how to deal with them
        update(currentCutoff);
        ClusterInformation data;
        const float dist = dMatrix->getSmallDist(); // Round to the third decimal place

        const float rndDist = Utils::ceilDist(dist, precision);
        if (previousDist <= 0.0000 && !Utils::isEqual(dist, previousDist)) {
            data.label = "0.00000";
            data.numberOfOtu = oldList.getNumBins();
        } else if (!Utils::isEqual(rndDist, rndPreviousDist)) {
            data.label = std::to_string(rndPreviousDist);
            data.numberOfOtu = oldList.getNumBins();
        }



        if(!data.label.empty()) {
            data.clusterBins = oldList.print(listFile);
            ListVector listVec(oldList);
            list->setPrintedLabels(false);
            clusterData->AddToData(data);
            if(rndPreviousDist > highestDistLabel) {
                highestDistLabel = rndPreviousDist;
                listVec.setLabel(std::to_string(highestDistLabel));
                clusterData->SetListVector(listVec, std::to_string(highestDistLabel)); // vec might be a shallow copy
            }
        }
        oldList = *list;
        previousDist = dist;
        rndPreviousDist = rndDist;
    }
    ClusterInformation data;
    if(previousDist <= 0.0000) {
        data.label = std::to_string(previousDist);
        data.numberOfOtu = oldList.getNumBins();
    }
    else if(rndPreviousDist<currentCutoff) {
        data.label = std::to_string(rndPreviousDist);
        data.numberOfOtu = oldList.getNumBins();
    }

    if(!data.label.empty()) {
        data.clusterBins = oldList.print(listFile);
        ListVector listVec(oldList);
        clusterData->AddToData(data);
        if(rndPreviousDist > highestDistLabel) {
            highestDistLabel = rndPreviousDist;
            listVec.setLabel(std::to_string(highestDistLabel));
            clusterData->SetListVector(listVec, std::to_string(highestDistLabel));
        }
    }
    return clusterData;
}

/***********************************************************************/
bool Cluster::update(double &cutOFF) {
    smallCol = dMatrix->getSmallestCell(smallRow);
    if (static_cast<int>(smallCol) == -1)
        return false;
    nColCells = dMatrix->seqVec[smallCol].size();
    nRowCells = dMatrix->seqVec[smallRow].size();
    std::vector<int> foundCol(nColCells, 0);

    bool changed = false;

    for (int i = static_cast<int>(nRowCells) - 1; i >= 0; i--) {
        //matrix indexes sorted from largest to smallest, so start at smallest index


        //if you are not the smallCell
        if (dMatrix->seqVec[smallRow][i].index != smallCol) {
            const long long search = dMatrix->seqVec[smallRow][i].index;

            bool merged = false;
            for (size_t j = 0; j < nColCells; j++) {
                //go through each distance the smallCol has looking for matching distance to find

                if (dMatrix->seqVec[smallCol][j].index != smallRow) {
                    //if you are not the smallest distance
                    if (dMatrix->seqVec[smallCol][j].index == search) {
                        //we found a distance for the merge
                        foundCol[j] = 1;
                        merged = true;
                        changed = updateDistance(dMatrix->seqVec[smallCol][j], dMatrix->seqVec[smallRow][i]);
                        dMatrix->updateCellCompliment(smallCol, j);
                        break;
                    } else if (dMatrix->seqVec[smallCol][j].index < search) {
                        //we don't have a distance for this cell
                        if (!Utils::isEqual(adjust, -1)) {
                            //adjust
                            merged = true;
                            const PDistCell value(search, adjust); //create a distance for the missing value
                            const int location = dMatrix->addCellSorted(smallCol, value);
                            changed = updateDistance(dMatrix->seqVec[smallCol][location], dMatrix->seqVec[smallRow][i]);
                            dMatrix->updateCellCompliment(smallCol, location);
                            nColCells++;
                            foundCol.push_back(0); //add a new-found column
                            //adjust value
                            for (int k = static_cast<int>(foundCol.size() - 1); k > location; k--) { foundCol[k] = foundCol[k - 1]; }
                            foundCol[location] = 1;
                        }
                        j += nColCells; //jump out of loop and remove cell below
                    }
                }
            }
            //if not merged it you need it for warning
            if ((!merged) && (method == "average" || method == "weighted")) {
                if (cutOFF > dMatrix->seqVec[smallRow][i].dist) { cutOFF = dMatrix->seqVec[smallRow][i].dist; }
            }
            if ((method == "nearest") && (!merged)) {
                //you are a row dist without a column dist, add you as a column dist
                const PDistCell value(search, dMatrix->seqVec[smallRow][i].dist); //create a distance for the missing value
                const int location = dMatrix->addCellSorted(smallCol, value);
                nColCells++;
                foundCol.push_back(0); //add a new-found column
                //adjust value
                for (int k = static_cast<int>(foundCol.size() - 1); k > location; k--) { foundCol[k] = foundCol[k - 1]; }
                foundCol[location] = 1;
            }
            dMatrix->rmCell(smallRow, i);
        }
    }
    clusterBins();
    clusterNames();

    if (method == "nearest") {
        for (int i = static_cast<int>(nColCells) - 1; i >= 0; i--) {
            //remove any unfound dists from merged column, need special case for nn, since unfound dists mean above the cutoff . keep smaller dist in col
            if (foundCol[i] == 0) {
                //not found
                if (dMatrix->seqVec[smallCol][i].index == smallRow) {
                    //you are the smallest distance
                    dMatrix->rmCell(smallCol, i);
                    break;
                }
            }
        }
    } else {
        for (int i = static_cast<int>(nColCells) - 1; i >= 0; i--) {
            //remove any unfound dists from merged column, need special case for nn, since unfound dists mean above the cutoff . keep smaller dist in col
            if (foundCol[i] == 0) {
                //not found
                if (!Utils::isEqual(adjust, -1)) {
                    //adjust
                    PDistCell value(static_cast<long long>(smallCol), adjust); //create a distance for the missing value
                    changed = updateDistance(dMatrix->seqVec[smallCol][i], value);
                    dMatrix->updateCellCompliment(smallCol, i);
                } else {
                    if (method == "average" || method == "weighted") {
                        if (dMatrix->seqVec[smallCol][i].index != smallRow) {
                            //if you are not hte the smallest distance
                            if (cutOFF > dMatrix->seqVec[smallCol][i].dist) {
                                cutOFF = dMatrix->seqVec[smallCol][i].dist;
                            }
                        }
                    }
                }
                dMatrix->rmCell(smallCol, i);
            }
        }
    }
    //dMatrix.print();
    return changed;
}

/***********************************************************************/
bool Cluster::updateMap() {
    //update location of seqs in smallRow since they move to smallCol now
    const std::string names = list->get(static_cast<int>(smallRow));
    std::string individual;
    const size_t binNameslength = names.size();
    if(names.empty())
        return false;

    for (size_t j = 0; j < binNameslength; j++) {
        if (names[j] == ',') {
            seq2Bin[individual] = static_cast<int>(smallCol);
            individual = "";
        } else { individual += names[j]; }
    }
    //get last name
    seq2Bin[individual] = static_cast<int>(smallCol);
    return true;
}

/***********************************************************************/



