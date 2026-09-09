#ifndef DOTAXONOMY_H
#define DOTAXONOMY_H

/*
 * phylotree.h
 *  
 *
 *  Created by Pat Schloss on 6/17/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include <map>
#include <set>
#include <utility>
#include <vector>
#include <string>
#include <fstream>
#include "Taxon.h"
/**************************************************************************************************/

struct TaxNode {
	std::vector<std::string> accessions;	//names of seqs in this branch of tree
	std::map<std::string, int> children;  //childs name to index in tree
	int parent, childNumber, level;
	std::string name, heirarchyID;
	
	TaxNode(std::string n) : parent(-1), level(0), name(std::move(n)) {		}
	TaxNode() = default;
};

/**************************************************************************************************/

class PhyloTree {

public:
	PhyloTree();
	PhyloTree(std::string);  //pass it a taxonomy file and it makes the tree
	PhyloTree(std::ifstream&, std::string);  //pass it a taxonomy file and it makes the train.tree
	~PhyloTree() = default;;
	int addSeqToTree(std::string, std::string);
    int addSeqToTree(const std::string &seqName, const std::vector<Taxon> &seqTax);
	void assignHeirarchyIDs(int);
	void printTreeNodes(std::string); //used by bayesian to save time
	std::vector<int> getGenusNodes();
	std::vector<int> getGenusTotals();	
	void setUp(std::string);  //used to create file needed for summary file if you use () constructor and add seqs manually instead of passing taxonomyfile
		
	TaxNode get(int i);				
	TaxNode get(std::string seqName);
    std::vector<TaxNode> getNodes(int); //returns std::vector of nodes at given level
	std::string getName(int i);			
	int getGenusIndex(std::string seqName);	
	std::string getFullTaxonomy(std::string);	 //pass a sequence name return taxonomy
    std::vector<std::string> getSeqs(std::string);	 //returns names of sequences in given taxonomy
	
	int getMaxLevel()		{	return maxLevel;	}
	int getNumSeqs()		{	return numSeqs;		}
	int getNumNodes()		{	return static_cast<int>(tree.size());	}
	
	bool ErrorCheck(std::vector<std::string>);
	
private:
	void print(std::ofstream&, std::vector<TaxNode>&); //used to create static reference taxonomy file
	void fillOutTree(int, std::vector<TaxNode>&); //used to create static reference taxonomy file
	void binUnclassified(std::string);
	
	std::vector<TaxNode> tree;
	std::vector<int> genusIndex; //holds the indexes in tree where the genus level taxonomies are stored
	std::vector<int> totals; //holds the numSeqs at each genus level taxonomy
	std::map<std::string, int> name2Taxonomy;  //std::maps name to index in tree
    std::map<std::string, int> name2GenusNodeIndex;
	std::set<int> uniqueTaxonomies;  //std::map of unique taxonomies
	std::map<int, int> leafNodes; //used to create static reference taxonomy file
	//void print(int, ofstream&);
	int numNodes;
	int numSeqs;
	int maxLevel;
	bool calcTotals;
};

/**************************************************************************************************/

#endif


