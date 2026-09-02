/*
 *  splitmatrix.cpp
 *  Mothur
 *
 *  Created by westcott on 5/19/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 */


#include "splitmatrix.h"
#include "../MothurDependencies/Utils.h"

/***********************************************************************/
SplitMatrix::SplitMatrix(std::string ffile, std::string name, std::string count, std::string tax, float c, float cu, int p, bool cl, std::string output, bool v){
	fastafile = ffile;
	namefile = name;
    countfile = count;
	taxFile = tax;
	cutoff = c;  //tax level cutoff
	distCutoff = cu; //for fasta method if you are creating distance matrix you need a cutoff for that
	processors = p;
    classic = cl;
	outputDir = output;
    usingVsearchToCLuster = v;
    
    splitClassify();
   
}

SplitMatrix::SplitMatrix(const std::vector<std::string> &fasta, const CountTableAdapter &countTable,
	const std::vector<TaxonomyData>& data, float c, float distCut, int threads, bool cl, bool v) {


}

/***********************************************************************/
void SplitMatrix::splitClassify(){
	cutoff = int(cutoff);
    std::map<std::string, std::string> temp; Utils::readTax(taxFile, temp, true);

    PhyloTree phylo;
    for (std::map<std::string, std::string>::iterator itTemp = temp.begin(); itTemp != temp.end();) {

        if (m->getControl_pressed()) { return; }

        phylo.addSeqToTree(itTemp->first, itTemp->second);
        temp.erase(itTemp++);
    }

	phylo.assignHeirarchyIDs(0);

	//make sure the cutoff is not greater than maxlevel
	if (cutoff > phylo.getMaxLevel()) { m->mothurOut("splitcutoff is greater than the longest taxonomy, using " + tostd::string(phylo.getMaxLevel())); m->mothurOutEndLine(); cutoff = phylo.getMaxLevel(); }

    std::vector<std::vector<std::string> > seqGroups; //seqFroups[0] -> std::vector of std::string containing names of seqs assigned to group 0
    std::vector<std::string> taxGroupNames;
	//for each node in tree
	for (int i = 0; i < phylo.getNumNodes(); i++) {

        if (m->getControl_pressed()) { return; }

		//is this node within the cutoff
		TaxNode taxon = phylo.get(i);

		if (taxon.level == cutoff) {//if yes, then create group containing this nodes sequences
			if (taxon.accessions.size() > 1) { //if this taxon just has one seq its a singleton
                std::vector<std::string> thisGroupsSeqs;
                for (int j = 0; j < taxon.accessions.size(); j++) {
                    thisGroupsSeqs.push_back(taxon.accessions[j]);
				}
                seqGroups.push_back(thisGroupsSeqs);
                taxGroupNames.push_back(taxon.name);
			}
		}
	}

    if (usingVsearchToCLuster)  { createFastaFilesFromTax(seqGroups, taxGroupNames);        }
    else                        {  createDistanceFilesFromTax(seqGroups, taxGroupNames);    }

}
/***********************************************************************/
int SplitMatrix::createDistanceFilesFromTax(std::vector<std::vector<std::string> >& seqGroups, std::vector<std::string> groupNames){
    try {
        
        int numGroups = seqGroups.size();
        std::string thisOutputDir = outputDir;
        if (outputDir == "") {  thisOutputDir += Utils::hasPath(fastafile);  }
        std::string nonSingletonsFile = thisOutputDir + Utils::getRootName(Utils::getSimpleName(fastafile)) + "nonsingleton.accnos";
        ofstream outNonSingleton; Utils::openOutputFile(nonSingletonsFile, outNonSingleton);
        
        ifstream inFASTA; Utils::openInputFile(fastafile, inFASTA);
        SequenceDB fullDB(inFASTA); inFASTA.close();
        
        if (m->getDebug()) { for (int i = 0; i < numGroups; i++) { m->mothurOut("[DEBUG]: Number of unique sequences for group " + groupNames[i] + " (" + tostd::string(i+1) + " of " + tostd::string(numGroups) + "): " + tostd::string(seqGroups[i].size()) + "\n\n"); } }
        
        //process each group
        for (int i = 0; i < numGroups; i++) {
            
            if (m->getControl_pressed()) { outNonSingleton.close(); Utils::mothurRemove(nonSingletonsFile); for (int i = 0; i < dists.size(); i++) { Utils::mothurRemove((dists[i].begin()->first)); Utils::mothurRemove((dists[i].begin()->second)); } dists.clear(); return 0; }
            
            unordered_set<std::string> thisGroupsNames = Utils::mothurConvert(seqGroups[i]);
            
            m->mothurOut("/******************************************/\n");
            m->mothurOut("Selecting sequences for group " + groupNames[i] + " (" + tostd::string(i+1) + " of " + tostd::string(numGroups) + ")\nNumber of unique sequences: " + tostd::string(seqGroups[i].size()) + "\n\n");
            
            std::string outName = "";
            if (namefile != "") {
                thisOutputDir = outputDir;
                if (outputDir == "") {  thisOutputDir += Utils::hasPath(namefile);  }
                outName = thisOutputDir + Utils::getRootName(Utils::getSimpleName(namefile)) + tostd::string(i) + ".name.temp";
            }
            
            if (countfile != "") {
                thisOutputDir = outputDir;
                if (outputDir == "") {  thisOutputDir += Utils::hasPath(countfile);  }
                outName = thisOutputDir + Utils::getRootName(Utils::getSimpleName(countfile)) + tostd::string(i) + ".count.temp";
            }
            
            pair<std::string,std::string> dupsFile(namefile, outName); std::string dupsFormat = "name";
            if (countfile != "") { dupsFile.first = countfile; dupsFormat = "count"; }
            
            Command* getCommand = new GetSeqsCommand(thisGroupsNames, nullstd::stringPair, nullstd::stringPair, dupsFile, dupsFormat);
            
            delete getCommand;
            
            StorageDatabase* thisDB;
            
            m->mothurOut("\nCalculating distances for group " + groupNames[i] + " (" + tostd::string(i+1) + " of " + tostd::string(numGroups) + "):\n");
            thisOutputDir = outputDir;
            if (outputDir == "") {  thisOutputDir += Utils::hasPath(fastafile);  }
            std::string outputFileRoot = thisOutputDir + Utils::getRootName(Utils::getSimpleName(fastafile)) + tostd::string(i) + ".";
            
            std::string outputformat = "column"; if (classic) { outputformat = "lt"; }
            
            Command* command;
            std::vector< std::vector< int > > kmerDB; std::vector< int > lengths;
            
            if (fullDB.sameLength()) {
                thisDB = new SequenceDB(fullDB, thisGroupsNames);
                command = new DistanceCommand(thisDB, outputFileRoot, distCutoff, outputformat, processors);
            }
            else {
                thisDB = new SequenceDB(fullDB, thisGroupsNames, 7, kmerDB, lengths);
                command = new PairwiseSeqsCommand(thisDB, kmerDB, lengths, outputFileRoot, distCutoff, outputformat, processors);
            }
            
            std::map<std::string, std::vector<std::string> > filenames = command->getOutputFiles();
            
            std::string thisDistanceFile = "";
            if (classic) { thisDistanceFile = filenames["phylip"][0]; }
            else { thisDistanceFile = filenames["column"][0]; }
            
            delete command;
            
            m->mothurOut("/******************************************/\n");
            
            if (!Utils::isBlank(thisDistanceFile)) {
                thisOutputDir = outputDir;
                if (outputDir == "") {  thisOutputDir += Utils::hasPath(fastafile);  }
                std::string outDist = thisOutputDir + Utils::getRootName(Utils::getSimpleName(fastafile)) + tostd::string(i) + ".disttemp";
                Utils::renameFile(thisDistanceFile, outDist);
                
                std::map<std::string, std::string> thisFilePair; thisFilePair[outDist] = outName;
                dists.push_back(thisFilePair);
                
                for (int j = 0; j < thisDB->getNumSeqs(); j++) { outNonSingleton << thisDB->getSeq(j).getName() << endl; }
            }else {
                Utils::mothurRemove(thisDistanceFile);
                Utils::mothurRemove(outName);
            }
            
            delete thisDB;
        }
        outNonSingleton.close();
        
        
        if (!Utils::isBlank(nonSingletonsFile)) { //there are non singletons, so remove them to find the singletons
            //get singletons
            if (namefile != "") {
                thisOutputDir = outputDir;
                if (outputDir == "") {  thisOutputDir += Utils::hasPath(namefile);  }
                singleton = thisOutputDir + Utils::getRootName(Utils::getSimpleName(namefile)) + "singletons.temp";
            }
            
            if (countfile != "") {
                thisOutputDir = outputDir;
                if (outputDir == "") {  thisOutputDir += Utils::hasPath(countfile);  }
                singleton = thisOutputDir + Utils::getRootName(Utils::getSimpleName(countfile)) + "singletons.temp";
            }
            
            pair<std::string,std::string> dupsFile(namefile, singleton); std::string dupsFormat = "name";
            if (countfile != "") { dupsFile.first = countfile; dupsFormat = "count"; }
            
            m->mothurOut("/******************************************/\n");
            m->mothurOut("Finding singletons (ignore 'Removing group' messages):\n\nRunning command: remove.seqs()\n");
            
            Command* removeCommand = new RemoveSeqsCommand(nonSingletonsFile, dupsFile, dupsFormat);
            
            delete removeCommand;
            
            m->mothurOut("/******************************************/\n");
            
        }else { //every seqs is a singleton
            if (namefile != "") { singleton = namefile; }
            else if (countfile != "") { singleton = countfile; }
        }
        
        if (Utils::isBlank(singleton)) {
            Utils::mothurRemove(singleton);
            singleton = "none";
        }
        
        Utils::mothurRemove(nonSingletonsFile);
        
        if (m->getControl_pressed())     {  for (int i = 0; i < dists.size(); i++) { Utils::mothurRemove((dists[i].begin()->first)); Utils::mothurRemove((dists[i].begin()->second)); } dists.clear(); }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "SplitMatrix", "createDistanceFilesFromTax");
        exit(1);
    }
}
/***********************************************************************/
int SplitMatrix::createFastaFilesFromTax(std::vector<std::vector<std::string> >& seqGroups, std::vector<std::string> groupNames){
    try {
        
        int numGroups = seqGroups.size();
        std::string thisOutputDir = outputDir;
        if (outputDir == "") {  thisOutputDir += Utils::hasPath(fastafile);  }
        singleton = "none";
          
        ifstream inFASTA; Utils::openInputFile(fastafile, inFASTA);
        SequenceDB fullDB(inFASTA); inFASTA.close();
        
        if (!fullDB.sameLength()) {
            m->mothurOut("[ERROR]: Cannot cluster using vsearch with unaligned sequences, please correct.\n\n"); m->setControl_pressed(true);
        }
        
        if (m->getDebug()) { for (int i = 0; i < numGroups; i++) { m->mothurOut("[DEBUG]: Number of unique sequences for group " + groupNames[i] + " (" + tostd::string(i+1) + " of " + tostd::string(numGroups) + "): " + tostd::string(seqGroups[i].size()) + "\n\n"); } }
        
        //process each group
        for (int i = 0; i < numGroups; i++) {
            
            if (m->getControl_pressed()) {  for (int i = 0; i < dists.size(); i++) { Utils::mothurRemove((dists[i].begin()->first)); Utils::mothurRemove((dists[i].begin()->second)); } dists.clear(); return 0; }
            
            unordered_set<std::string> thisGroupsNames = Utils::mothurConvert(seqGroups[i]);
            
            m->mothurOut("/******************************************/\n");
            m->mothurOut("Selecting sequences for group " + groupNames[i] + " (" + tostd::string(i+1) + " of " + tostd::string(numGroups) + ")\nNumber of unique sequences: " + tostd::string(seqGroups[i].size()) + "\n\n");
            
            std::string outName = "";
            if (namefile != "") {
                thisOutputDir = outputDir;
                if (outputDir == "") {  thisOutputDir += Utils::hasPath(namefile);  }
                outName = thisOutputDir + Utils::getRootName(Utils::getSimpleName(namefile)) + tostd::string(i) + ".name.temp";
            }
            
            if (countfile != "") {
                thisOutputDir = outputDir;
                if (outputDir == "") {  thisOutputDir += Utils::hasPath(countfile);  }
                outName = thisOutputDir + Utils::getRootName(Utils::getSimpleName(countfile)) + tostd::string(i) + ".count.temp";
            }
            
            pair<std::string, std::string> dupsFile(namefile, outName); std::string dupsFormat = "name";
            if (countfile != "") { dupsFile.first = countfile; dupsFormat = "count"; }
            
            Command* getCommand = new GetSeqsCommand(thisGroupsNames, nullstd::stringPair, nullstd::stringPair, dupsFile, dupsFormat);
            
            delete getCommand;
            
            StorageDatabase* thisDB;
            
            thisOutputDir = outputDir;
            if (outputDir == "") {  thisOutputDir += Utils::hasPath(fastafile);  }
            std::string outFasta = thisOutputDir + Utils::getRootName(Utils::getSimpleName(fastafile)) + tostd::string(i) + ".fastatemp";
            thisDB = new SequenceDB(fullDB, thisGroupsNames);
            thisDB->print(outFasta);
            delete thisDB;
            
            std::map<std::string, std::string> thisFilePair; thisFilePair[outFasta] = outName;
            dists.push_back(thisFilePair);
        }
        
        if (m->getControl_pressed())     {  for (int i = 0; i < dists.size(); i++) { Utils::mothurRemove((dists[i].begin()->first)); Utils::mothurRemove((dists[i].begin()->second)); } dists.clear(); }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "SplitMatrix", "createDistanceFilesFromTax");
        exit(1);
    }
}
/********************************************************************************************************************/
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
std::vector< std::map< std::string, std::string> > SplitMatrix::getDistanceFiles(){
	try {	
		
		sort(dists.begin(), dists.end(), compareFileSizes);
		
		return dists;
	}
	catch(exception& e) {
		m->errorOut(e, "SplitMatrix", "getDistanceFiles");
		exit(1);
	}
}
//********************************************************************************************************************/


