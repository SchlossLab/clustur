#ifndef SPLITMATRIX_H
#define SPLITMATRIX_H
/*
 *  splitmatrix.h
 *  Mothur
 *
 *  Created by westcott on 5/19/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

/******************************************************/
#include <vector>
#include <string>
#include <map>

#include "OptidataComponents.h"
#include "../DataStructures/TaxonomyData.h"

#include "../Adapters/CountTableAdapter.h"
#include "../MothurDependencies/PairwiseDistanceCalculator.h"

class SplitMatrix  {
	
public:

 //    SplitMatrix(std::string ffile, std::string name, std::string count, std::string tax,
 //    	float c, float cu, int p, bool cl, std::string output, bool v); //fastafile, namefile, countfile, taxFile, taxcutoff, cutoff, processors, classic, outputDir, usingVsearchToCLuster
	// SplitMatrix(const std::vector<TaxonomyData>& data, const FastaDatabase& fastaDatabase,
	// PairwiseDistanceCalculator* calculator, float c, int taxLevel);


	static std::vector<OptiDataComponent> splitClassify(const std::vector<TaxonomyData>& data, const FastaDatabase& fastaDatabase,
	PairwiseDistanceCalculator* calculator, double cutoff, int taxLevel);


    // SplitMatrix(const std::vector<std::string> &fasta, const CountTableAdapter &countTable,
    //             const std::vector<TaxonomyData>& data, float c, float distCut, int threads, bool cl, bool v = false);
    ~SplitMatrix() = default;
    
    // std::vector< std::map<std::string, std::string> > getDistanceFiles();  //returns map of distance files -> namefile sorted by distance file size

    // std::string getSingletonNames() { return singleton; } //returns namesfile or countfile containing singletons
    //long long getNumSingleton() { return numSingleton; } //returns namesfile containing singletons
	
private:
    // std::string distFile, namefile, singleton,  taxFile, fastafile, outputDir, countfile;
    // std::vector< std::map< std::string, std::string> > dists;
    // float cutoff, distCutoff;
    // bool classic, usingVsearchToCLuster;
    // int processors;
	// PairwiseDistanceCalculator* calculator;
	// std::vector<OptiData*> optidataMatrices;

	static std::vector<OptiDataComponent> createDistanceFilesFromTax(const FastaDatabase& fastaData, PairwiseDistanceCalculator* calculator,
		std::vector<std::vector<std::string> >& seqGroups, const std::vector<std::string> &groupNames, double cutoff);
    // int createDistanceFilesFromTax(std::vector<std::vector<std::string> >&, std::vector<std::string>);
    // int createFastaFilesFromTax(std::vector<std::vector<std::string> >&, std::vector<std::string>);
};

/******************************************************/

#endif

