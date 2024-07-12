//
// Created by Gregory Johnson on 6/17/24.
//

#include "MothurDependencies/RAbundVector.h"
/*
 *  rabundvector.cpp
 *
 *
 *  Created by Pat Schloss on 8/8/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

/***********************************************************************/

RAbundVector::RAbundVector() : DataVector(), maxRank(0), numBins(0), numSeqs(0) {}

/***********************************************************************/

RAbundVector::RAbundVector(int n) : DataVector(), data(n,0) , maxRank(0), numBins(0), numSeqs(0) {}


/***********************************************************************/

RAbundVector::~RAbundVector() {

}

/***********************************************************************/

void RAbundVector::set(int binNumber, int newBinSize){
	int oldBinSize = data[binNumber];
	data[binNumber] = newBinSize;

	if(oldBinSize == 0)			{	numBins++;				}
	if(newBinSize == 0)			{	numBins--;				}
	if(newBinSize > maxRank)	{	maxRank = newBinSize;	}

	numSeqs += (newBinSize - oldBinSize);
}

/***********************************************************************/

int RAbundVector::get(const int index) const {
	return data[index];

}
/***********************************************************************/

void RAbundVector::clear(){
	numBins = 0;
	maxRank = 0;
	numSeqs = 0;
	data.clear();

}
/***********************************************************************/

void RAbundVector::push_back(const int binSize){
		data.push_back(binSize);
		numBins++;

		if(binSize > maxRank){
			maxRank = binSize;
		}

		numSeqs += binSize;
}
/***********************************************************************/

int RAbundVector::remove(const int bin){
    const int abund = data[bin];
    data.erase(data.begin()+bin);
    numBins--;

    if(abund == maxRank){
        const auto it = max_element(data.begin(), data.end());
        maxRank = *it;
    }

    numSeqs -= abund;

    return abund;

}
/***********************************************************************/

void RAbundVector::resize(int size){

	data.resize(size);
}

/***********************************************************************/

int RAbundVector::size(){
	return data.size();
}

/***********************************************************************/

bool RAbundVector::quicksort(){
	std::sort(data.rbegin(), data.rend());
	return true;
}


/***********************************************************************/
int RAbundVector::getNumBins(){
	return numBins;
}

/***********************************************************************/

int RAbundVector::getNumSeqs(){
	return numSeqs;
}

/***********************************************************************/

int RAbundVector::getMaxRank(){
	return maxRank;
}

/***********************************************************************/

// RAbundFloatVector RAbundVector::getRAbundFloatVector(){
//     RAbundFloatVector rav; rav.setLabel(label);
//     for(int i=0;i<data.size();i++){ rav.push_back(0.0 + data[i]);  }
//     return rav;
// }

// /***********************************************************************/
//
// SAbundVector RAbundVector::getSAbundVector() {
// 	try {
// 		SAbundVector sav(maxRank+1);
//
// 		for(int i=0;i<data.size();i++){
// 			int abund = data[i];
// 			sav.set(abund, sav.get(abund) + 1);
// 		}
// 		sav.set(0, 0);
// 		sav.setLabel(label);
// 		return sav;
// 	}
// 	catch(exception& e) {
// 		m->errorOut(e, "RAbundVector", "getSAbundVector");
// 		exit(1);
// 	}
// }

/***********************************************************************/

// OrderVector RAbundVector::getOrderVector(map<std::string,int>* nameMap = nullptr) {
// 	try {
//
//         std::vector<int> ovData;
// 		for(int i=0;i<data.size();i++){
// 			for(int j=0;j<data[i];j++){
// 				ovData.push_back(i);
// 			}
// 		}
//
// 		util.mothurRandomShuffle(ovData);
//
//         OrderVector ov;
//         for(int i=0;i<ovData.size();i++){ ov.push_back(ovData[i]); }
// 		ov.setLabel(label);
// 		ov.getNumBins();
//
// 		return ov;
// 	}
// 	catch(exception& e) {
// 		m->errorOut(e, "RAbundVector", "getOrderVector");
// 		exit(1);
// 	}
// }

/***********************************************************************/
