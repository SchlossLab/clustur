//
// Created by Gregory Johnson on 6/17/24.
//

#include "DataStructures/RAbundVector.h"
/*
 *  rabundvector.cpp
 *
 *
 *  Created by Pat Schloss on 8/8/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

/***********************************************************************/

RAbundVector::RAbundVector() : maxRank(0), numBins(0), numSeqs(0) {}

/***********************************************************************/

RAbundVector::RAbundVector(const int n) : data(n,0) , maxRank(0), numBins(0), numSeqs(0) {}


/***********************************************************************/

RAbundVector::~RAbundVector() {

}

/***********************************************************************/

void RAbundVector::set(const int binNumber, const int newBinSize){
	const int oldBinSize = data[binNumber];
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

void RAbundVector::resize(const int size){

	data.resize(size);
}

/***********************************************************************/

int RAbundVector::size() const {
	return static_cast<int>(data.size());
}

/***********************************************************************/

bool RAbundVector::quicksort(){
	std::sort(data.rbegin(), data.rend());
	return true;
}


/***********************************************************************/
int RAbundVector::getNumBins() const {
	return numBins;
}

/***********************************************************************/

int RAbundVector::getNumSeqs() const {
	return numSeqs;
}

/***********************************************************************/

int RAbundVector::getMaxRank() const {
	return maxRank;
}

/***********************************************************************/
void RAbundVector::print() const {
	const std::string path = "/Users/grejoh/Documents/OptiClusterPackage/clustur/output_rabund.txt";
	std::ofstream output(path);
	output << label;
	output << '\t' << numBins;

	std::vector<int> hold = data;
	std::sort(hold.rbegin(), hold.rend());

	for(int i=0;i<numBins;i++){	output  << '\t' << hold[i];	}

	output << std::endl;

	output.close();
}
