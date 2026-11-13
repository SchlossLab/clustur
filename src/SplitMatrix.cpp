/*
 *  splitmatrix.cpp
 *  Mothur
 *
 *  Created by westcott on 5/19/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 */


#include "DataStructures/FastaDatabase.h"
#include "DataStructures/OptiData.h"
#include "DataStructures/OptiDataComponentFactory.h"
#include "DataStructures/PhyloTree.h"
#include "DataStructures/SplitMatrix.h"
#include "MothurDependencies/OneGapPairwiseDistance.h"
#include "MothurDependencies/PairwiseDistanceCalculator.h"
#include "MothurDependencies/Utils.h"

/***********************************************************************/
// SplitMatrix::SplitMatrix(std::string ffile, std::string name, std::string count, std::string tax, float c, float cu, int p, bool cl, std::string output, bool v){
// 	fastafile = ffile;
// 	namefile = name;
//     countfile = count;
// 	taxFile = tax;
// 	cutoff = c;  //tax level cutoff
// 	distCutoff = cu; //for fasta method if you are creating distance matrix you need a cutoff for that
// 	processors = p;
//     classic = cl;
// 	outputDir = output;
//     usingVsearchToCLuster = v;
//
//     //splitClassify(TODO);
//
// }

// SplitMatrix::SplitMatrix(const std::vector<TaxonomyData>& data, const FastaDatabase& fastaDatabase,
// 	PairwiseDistanceCalculator* calculator, const float c,
// 	const int taxLevel):
// calculator(calculator){
// 	// fastafile = ffile;
// 	// namefile = name;
// 	// countfile = count;
// 	// taxFile = tax;
// 	cutoff = c;  //tax level cutoff
// 	// distCutoff = cu; //for fasta method if you are creating distance matrix you need a cutoff for that
// 	// processors = p;
// 	// classic = cl;
// 	// outputDir = output;
// 	// usingVsearchToCLuster = v;
// 	splitClassify(data, fastaDatabase);
// }

/***********************************************************************/
// void SplitMatrix::splitClassify(){
// 	cutoff = int(cutoff);
//     std::map<std::string, std::string> temp; Utils::readTax(taxFile, temp, true);
//
//     PhyloTree phylo;
//     for (std::map<std::string, std::string>::iterator itTemp = temp.begin(); itTemp != temp.end();) {
//
//         if (m->getControl_pressed()) { return; }
//
//         phylo.addSeqToTree(itTemp->first, itTemp->second);
//         temp.erase(itTemp++);
//     }
//
// 	phylo.assignHeirarchyIDs(0);
//
// 	//make sure the cutoff is not greater than maxlevel
// 	if (cutoff > phylo.getMaxLevel()) {mothurOut("splitcutoff is greater than the longest taxonomy, using " + tostd::string(phylo.getMaxLevel())); m->mothurOutEndLine(); cutoff = phylo.getMaxLevel(); }
//
//     std::vector<std::vector<std::string> > seqGroups; //seqFroups[0] -> std::vector of std::string containing names of seqs assigned to group 0
//     std::vector<std::string> taxGroupNames;
// 	//for each node in tree
// 	for (int i = 0; i < phylo.getNumNodes(); i++) {
//
//         if (m->getControl_pressed()) { return; }
//
// 		//is this node within the cutoff
// 		TaxNode taxon = phylo.get(i);
//
// 		if (taxon.level == cutoff) {//if yes, then create group containing this nodes sequences
// 			if (taxon.accessions.size() > 1) { //if this taxon just has one seq its a singleton
//                 std::vector<std::string> thisGroupsSeqs;
//                 for (int j = 0; j < taxon.accessions.size(); j++) {
//                     thisGroupsSeqs.push_back(taxon.accessions[j]);
// 				}
//                 seqGroups.push_back(thisGroupsSeqs);
//                 taxGroupNames.push_back(taxon.name);
// 			}
// 		}
// 	}
//
//     if (usingVsearchToCLuster)  { createFastaFilesFromTax(seqGroups, taxGroupNames);        }
//     else                        {  createDistanceFilesFromTax(seqGroups, taxGroupNames);    }
//
// }

std::vector<OptiDataComponent>  SplitMatrix::splitClassify(const std::vector<TaxonomyData>& data, const FastaDatabase& fastaDatabase,
	const PairwiseDistanceCalculator* calculator, const double cutoff, const int taxLevel, const int numberOfThreads) {
	PhyloTree phylo;
	for (const auto&[name, taxonomy] : data) {
		phylo.addSeqToTree(name, taxonomy);
	}
	phylo.assignHeirarchyIDs(0);
	if (taxLevel > phylo.getMaxLevel()) {
		// Rcpp::message("splitcutoff is greater than the longest taxonomy, using " +
		// std::to_string(phylo.getMaxLevel())); taxLevel = phylo.getMaxLevel();
	}
	std::vector<std::vector<std::string> > seqGroups; //seqFroups[0] -> vector of string containing names of seqs assigned to group 0
	std::vector<std::string> taxGroupNames;
	//for each node in tree
	for (int i = 0; i < phylo.getNumNodes(); i++) {

		//is this node within the cutoff

		if (TaxNode taxon = phylo.get(i); taxon.level == taxLevel) {//if yes, then create group containing this nodes sequences
			if (taxon.accessions.size() > 1) { //if this taxon just has one seq its a singleton
				std::vector<std::string> thisGroupsSeqs;
				for (const auto & accession : taxon.accessions) {
					thisGroupsSeqs.push_back(accession);
				}
				seqGroups.push_back(thisGroupsSeqs);
				taxGroupNames.push_back(taxon.name);
			}
		}
	}
	// if (usingVsearchToCLuster)  { createFastaFilesFromTax(seqGroups, taxGroupNames);        }
	// else                        {  createDistanceFilesFromTax(seqGroups, taxGroupNames);    }
	return createDistanceFilesFromTax(fastaDatabase, calculator, seqGroups, taxGroupNames, cutoff, numberOfThreads);
}

/***********************************************************************/
std::vector<OptiDataComponent>  SplitMatrix::createDistanceFilesFromTax(const FastaDatabase& fastaData, const PairwiseDistanceCalculator* calculator,
	std::vector<std::vector<std::string> >& seqGroups, const std::vector<std::string> &groupNames, const double cutoff,
	const int numberOfThreads) {

    const size_t numGroups = seqGroups.size();
    //process each group
	std::vector<OptiDataComponent> result(numGroups);
    for (size_t i = 0; i < numGroups; i++) {
		std::sort(seqGroups[i].begin(), seqGroups[i].end());
    	result[i] = OptiDataComponentFactory::CreateOptiDataComponent(fastaData, calculator, seqGroups[i],
    		cutoff, numberOfThreads);
    }
    return result;
}
// /***********************************************************************/
// int SplitMatrix::createFastaFilesFromTax(std::vector<std::vector<std::string> >& seqGroups, std::vector<std::string> groupNames){
//     try {
//
//         int numGroups = seqGroups.size();
//         std::string thisOutputDir = outputDir;
//         if (outputDir == "") {  thisOutputDir += Utils::hasPath(fastafile);  }
//         singleton = "none";
//
//         ifstream inFASTA; Utils::openInputFile(fastafile, inFASTA);
//         SequenceDB fullDB(inFASTA); inFASTA.close();
//
//         if (!fullDB.sameLength()) {
//             m->mothurOut("[ERROR]: Cannot cluster using vsearch with unaligned sequences, please correct.\n\n"); m->setControl_pressed(true);
//         }
//
//         if (m->getDebug()) { for (int i = 0; i < numGroups; i++) { m->mothurOut("[DEBUG]: Number of unique sequences for group " + groupNames[i] + " (" + tostd::string(i+1) + " of " + tostd::string(numGroups) + "): " + tostd::string(seqGroups[i].size()) + "\n\n"); } }
//
//         //process each group
//         for (int i = 0; i < numGroups; i++) {
//
//             if (m->getControl_pressed()) {  for (int i = 0; i < dists.size(); i++) { Utils::mothurRemove((dists[i].begin()->first)); Utils::mothurRemove((dists[i].begin()->second)); } dists.clear(); return 0; }
//
//             unordered_set<std::string> thisGroupsNames = Utils::mothurConvert(seqGroups[i]);
//
//             m->mothurOut("/******************************************/\n");
//             m->mothurOut("Selecting sequences for group " + groupNames[i] + " (" + tostd::string(i+1) + " of " + tostd::string(numGroups) + ")\nNumber of unique sequences: " + tostd::string(seqGroups[i].size()) + "\n\n");
//
//             std::string outName = "";
//             if (namefile != "") {
//                 thisOutputDir = outputDir;
//                 if (outputDir == "") {  thisOutputDir += Utils::hasPath(namefile);  }
//                 outName = thisOutputDir + Utils::getRootName(Utils::getSimpleName(namefile)) + tostd::string(i) + ".name.temp";
//             }
//
//             if (countfile != "") {
//                 thisOutputDir = outputDir;
//                 if (outputDir == "") {  thisOutputDir += Utils::hasPath(countfile);  }
//                 outName = thisOutputDir + Utils::getRootName(Utils::getSimpleName(countfile)) + tostd::string(i) + ".count.temp";
//             }
//
//             pair<std::string, std::string> dupsFile(namefile, outName); std::string dupsFormat = "name";
//             if (countfile != "") { dupsFile.first = countfile; dupsFormat = "count"; }
//
//             Command* getCommand = new GetSeqsCommand(thisGroupsNames, nullstd::stringPair, nullstd::stringPair, dupsFile, dupsFormat);
//
//             delete getCommand;
//
//             StorageDatabase* thisDB;
//
//             thisOutputDir = outputDir;
//             if (outputDir == "") {  thisOutputDir += Utils::hasPath(fastafile);  }
//             std::string outFasta = thisOutputDir + Utils::getRootName(Utils::getSimpleName(fastafile)) + tostd::string(i) + ".fastatemp";
//             thisDB = new SequenceDB(fullDB, thisGroupsNames);
//             thisDB->print(outFasta);
//             delete thisDB;
//
//             std::map<std::string, std::string> thisFilePair; thisFilePair[outFasta] = outName;
//             dists.push_back(thisFilePair);
//         }
//
//         if (m->getControl_pressed())     {  for (int i = 0; i < dists.size(); i++) { Utils::mothurRemove((dists[i].begin()->first)); Utils::mothurRemove((dists[i].begin()->second)); } dists.clear(); }
//
//         return 0;
//     }
//     catch(exception& e) {
//         m->errorOut(e, "SplitMatrix", "createDistanceFilesFromTax");
//         exit(1);
//     }
// }
// /********************************************************************************************************************/
//sorts biggest to smallest
inline bool compareFileSizes(std::map<std::string, std::string> left, std::map<std::string, std::string> right){

	FILE * pFile;
	long leftsize = 0;

	//get num bytes in file
	std::string filename = left.begin()->first;
	pFile = fopen (filename.c_str(),"rb");
	std::string error = "Error opening " + filename;
	if (pFile==nullptr) perror (error.c_str());
	else{
		fseek (pFile, 0, SEEK_END);
		leftsize=ftell (pFile);
		fclose (pFile);
	}

	FILE * pFile2;
	long rightsize = 0;

	//get num bytes in file
	filename = right.begin()->first;
	pFile2 = fopen (filename.c_str(),"rb");
	error = "Error opening " + filename;
	if (pFile2==nullptr) perror (error.c_str());
	else{
		fseek (pFile2, 0, SEEK_END);
		rightsize=ftell (pFile2);
		fclose (pFile2);
	}

	return (leftsize > rightsize);
}
/***********************************************************************/
//returns std::map of distance files -> namefile sorted by distance file size
// std::vector< std::map< std::string, std::string> > SplitMatrix::getDistanceFiles(){
// 	sort(dists.begin(), dists.end(), compareFileSizes);
// 	return dists;
// }
