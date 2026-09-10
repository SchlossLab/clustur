/*
 *  doTaxonomy.cpp
 *
 *
 *  Created by Pat Schloss on 6/17/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "DataStructures/phylotree.h"
#include "MothurDependencies/Utils.h"

/**************************************************************************************************/

PhyloTree::PhyloTree(){
		numNodes = 1;
		numSeqs = 0;
		tree.emplace_back(TaxNode("Root"));
		tree[0].heirarchyID = "0";
        tree[0].level = 0;
		maxLevel = 0;
		calcTotals = true;
		addSeqToTree("unknown", "unknown;");
}
/**************************************************************************************************/

std::vector<std::string> PhyloTree::getSeqs(std::string seqTaxonomy){

    std::string taxCopy = seqTaxonomy;
    std::vector<std::string> names;
    std::map<std::string, int>::iterator childPointer;

	int currentNode = 0;

    Utils::removeConfidences(seqTaxonomy);
    std::vector<std::string> taxons; Utils::splitAtChar(seqTaxonomy, taxons, ';');

    for(std::string taxon : taxons) {

	//	if (m->getControl_pressed()) { return names; }

    //    if (m->getDebug()) { m->mothurOut(taxon +'\n'); }

	//	if (taxon == "") {  m->mothurOut(taxCopy + " has an error in the taxonomy.  This may be due to a ;;\n"); break;  }

		childPointer = tree[currentNode].children.find(taxon);

		if(childPointer != tree[currentNode].children.end()){	//if the node already exists, move on
			currentNode = childPointer->second;
		}
		else{											//otherwise, error this taxonomy is not in tree
	//		m->mothurOut("[ERROR]: " + taxCopy + " is not in taxonomy tree, please correct.\n"); m->setControl_pressed(true); return names;
		}
	}

    //return names in this taxonomy
    names = tree[currentNode].accessions;

    return names;
}
/**************************************************************************************************/

int PhyloTree::addSeqToTree(const std::string &seqName, const std::vector<Taxon> &seqTax){

    numSeqs++;

    std::map<std::string, int>::iterator childPointer;

    int currentNode = 0;
    int level = 0;

    tree[0].accessions.push_back(seqName);

    for (int i = 0; i < seqTax.size(); i++) {

        level++;
        std::string taxon = seqTax[i].name;


        childPointer = tree[currentNode].children.find(taxon);

        if(childPointer != tree[currentNode].children.end()){    //if the node already exists, move on
            currentNode = childPointer->second;
            tree[currentNode].accessions.push_back(seqName);
            name2Taxonomy[seqName] = currentNode;
        }
        else{                                            //otherwise, create it
            tree.emplace_back(taxon);

            tree[currentNode].children[taxon] = numNodes;
            tree[numNodes].level = level;
            tree[numNodes].parent = currentNode;
            currentNode = numNodes;
            numNodes++;

            tree[currentNode].accessions.push_back(seqName);
            name2Taxonomy[seqName] = currentNode;
        }

        if (i == (seqTax.size()-1)) {   uniqueTaxonomies.insert(currentNode);    } //last level
    }

    //save maxLevel for binning the unclassified seqs
    if (level > maxLevel) { maxLevel = level; }

    return level;
}
/**************************************************************************************************/

int PhyloTree::addSeqToTree(std::string seqName, std::string seqTaxonomy){

	numSeqs++;

	std::map<std::string, int>::iterator childPointer;

	tree[0].accessions.push_back(seqName);
	Utils::removeConfidences(seqTaxonomy);
    std::vector<std::string> taxons; Utils::splitAtChar(seqTaxonomy, taxons, ';');

    int level = 0;
    int currentNode = 0;

    for(std::string taxon : taxons) {

        level++;

	//	if (m->getControl_pressed()) { return 0; }

      //  if (m->getDebug()) { m->mothurOut(seqName +'\t' + taxon +'\n'); }

		//if (taxon == "") {  m->mothurOut(seqName + " has an error in the taxonomy.  This may be due to a ;;\n");  if (currentNode != 0) {  uniqueTaxonomies.insert(currentNode); } break;  }

		childPointer = tree[currentNode].children.find(taxon);

		if(childPointer != tree[currentNode].children.end()){	//if the node already exists, move on
			currentNode = childPointer->second;
			tree[currentNode].accessions.push_back(seqName);
			name2Taxonomy[seqName] = currentNode;
		}
		else{											//otherwise, create it
			tree.push_back(TaxNode(taxon));

			tree[currentNode].children[taxon] = numNodes;
            tree[numNodes].level = level;
			tree[numNodes].parent = currentNode;
            currentNode = numNodes;
            numNodes++;

			tree[currentNode].accessions.push_back(seqName);
			name2Taxonomy[seqName] = currentNode;
		}
	}

    uniqueTaxonomies.insert(currentNode);

    //save maxLevel for binning the unclassified seqs
    if (level > maxLevel) { maxLevel = level; }

	return level;
}
/**************************************************************************************************/
std::vector<int> PhyloTree::getGenusNodes()	{

	genusIndex.clear();
	//generate genusIndexes
	std::set<int>::iterator it2;
    std::map<int, int> temp;
	for (it2=uniqueTaxonomies.begin(); it2!=uniqueTaxonomies.end(); it2++) {  genusIndex.push_back(*it2); 	temp[*it2] = genusIndex.size()-1; }

    for (std::map<std::string, int>::iterator itName = name2Taxonomy.begin(); itName != name2Taxonomy.end(); itName++) {
        std::map<int, int>::iterator itTemp = temp.find(itName->second);
        if (itTemp != temp.end()) { name2GenusNodeIndex[itName->first] = itTemp->second; }
      //  else {  m->mothurOut("[ERROR]: trouble making name2GenusNodeIndex, aborting.\n"); m->setControl_pressed(true); }
    }

	return genusIndex;
}
/**************************************************************************************************/
std::vector<int> PhyloTree::getGenusTotals()	{

		if (calcTotals) {
			totals.clear();
			//reset counts because we are on a new word
			for (int j = 0; j < genusIndex.size(); j++) {
				totals.push_back(tree[genusIndex[j]].accessions.size());
			}
			return totals;
		}else{
			return totals;
		}

}
/**************************************************************************************************/
std::vector<TaxNode> PhyloTree::getNodes(int level)	{

    std::vector<TaxNode> thisLevelsNodes;

    //check level is not above max
    if (level > maxLevel) {
      //  m->mothurOut("[ERROR]: Tax level not present in file. Cannot continue.\n");
      //  m->setControl_pressed(true); return thisLevelsNodes;
    }

    for (size_t i = 0; i < tree.size(); i++) {
        if (tree[i].level == level) {
            thisLevelsNodes.push_back(tree[i]);
        }
    }

    return thisLevelsNodes;
}
/**************************************************************************************************/

void PhyloTree::assignHeirarchyIDs(int index){
	std::map<std::string,int>::iterator it;
	int counter = 1;

	for(it=tree[index].children.begin();it!=tree[index].children.end();it++){

    //    if (m->getDebug()) { m->mothurOut(tostd::string(index) +'\t' + tree[it->second].name +'\n'); }

		tree[it->second].heirarchyID = tree[index].heirarchyID + '.' + std::to_string(counter);
		counter++;
		//tree[it->second].level = tree[index].level + 1;

		//save maxLevel for binning the unclassified seqs
		if (tree[it->second].level > maxLevel) { maxLevel = tree[it->second].level; }

		assignHeirarchyIDs(it->second);
	}

}
/**************************************************************************************************/
void PhyloTree::setUp(std::string tfile){
	std::string taxFileNameTest = tfile.substr(0,tfile.find_last_of(".")+1) + "tree.sum";
    binUnclassified(taxFileNameTest);
}
/**************************************************************************************************/
void PhyloTree::binUnclassified(std::string file){

	std::ofstream out;
	//Utils::openOutputFile(file, out);

	std::map<std::string, int>::iterator itBin;
	std::map<std::string, int>::iterator childPointer;

	std::vector<TaxNode> copy = tree;

	//fill out tree
	fillOutTree(0, copy);

	//get leaf nodes that may need extension
	for (int i = 0; i < copy.size(); i++) {

		if (copy[i].children.size() == 0) {
			leafNodes[i] = i;
		}
	}

  //  if (m->getDebug()) { m->mothurOut("maxLevel = " + tostd::string(maxLevel) +'\n'); }

	int copyNodes = copy.size();

	//go through the seqs and if a sequence finest taxon is not the same level as the most finely defined taxon then classify it as unclassified where necessary
	std::map<int, int>::iterator itLeaf;
	for (itLeaf = leafNodes.begin(); itLeaf != leafNodes.end(); itLeaf++) {

	//	if (m->getControl_pressed()) {  out.close(); break;  }

		int level = copy[itLeaf->second].level;
		int currentNode = itLeaf->second;

     //   if (m->getDebug()) { m->mothurOut(copy[currentNode].name +'\n'); }

		//this sequence is unclassified at some levels
		while(level < maxLevel){

			level++;
         //   if (m->getDebug()) { m->mothurOut("level = " + tostd::string(level) +'\n'); }

			std::string taxon = "unclassified";

			//does the parent have a child names 'unclassified'?
			childPointer = copy[currentNode].children.find(taxon);

			if(childPointer != copy[currentNode].children.end()){	//if the node already exists, move on
				currentNode = childPointer->second; //currentNode becomes 'unclassified'
			}
			else{											//otherwise, create it
				copy.push_back(TaxNode(taxon));

				copy[currentNode].children[taxon] = copyNodes;
				copy[copyNodes].parent = currentNode;
				copy[copyNodes].level = copy[currentNode].level + 1;
                currentNode = copyNodes;
                copyNodes++;
			}
		}
	}

	// if (!m->getControl_pressed()) {
	// 	//print copy tree
	// 	print(out, copy);
	// }
}
/**************************************************************************************************/
void PhyloTree::fillOutTree(int index, std::vector<TaxNode>& copy) {

	std::map<std::string,int>::iterator it;

	it = copy[index].children.find("unclassified");
	if (it == copy[index].children.end()) { //no unclassified at this level
		std::string taxon = "unclassified";
		copy.push_back(TaxNode(taxon));
		copy[index].children[taxon] = copy.size()-1;
		copy[copy.size()-1].parent = index;
		copy[copy.size()-1].level = copy[index].level + 1;
	}

	if (tree[index].level < maxLevel) {
		for(it=tree[index].children.begin();it!=tree[index].children.end();it++){ //check your children
			fillOutTree(it->second, copy);
		}
	}
}
/**************************************************************************************************/
std::string PhyloTree::getFullTaxonomy(std::string seqName) {
	std::string tax = "";

	int currentNode = name2Taxonomy[seqName];

	while (tree[currentNode].parent != -1) {
		tax = tree[currentNode].name + ";" + tax;
		currentNode = tree[currentNode].parent;
	}

	return tax;
}
/**************************************************************************************************/

void PhyloTree::print(std::ofstream& out, std::vector<TaxNode>& copy){


	// //output mothur version
	// out << "#" << current->getVersion() << endl;
 //
	// out << copy.size() << endl;
 //
	// out << maxLevel << endl;
 //
	// for (int i = 0; i < copy.size(); i++) {
 //
	// 	out << copy[i].level << '\t' << copy[i].children.size() << endl;
 //        out << copy[i].name << endl;
 //
	// 	std::map<std::string,int>::iterator it;
	// 	for(it=copy[i].children.begin();it!=copy[i].children.end();it++){
	// 		out << it->second << '\t' << it->first << endl;
	// 	}
	// 	out << endl;
	// }
 //
	// out.close();
}
/**************************************************************************************************/
void PhyloTree::printTreeNodes(std::string treefilename) {

    // ofstream outTree;
    // Utils::openOutputFile(treefilename, outTree);
    //
    // //output mothur version
    // outTree << "#" << current->getVersion() << endl;
    //
    // //print treenodes
    // outTree << tree.size() << endl;
    //
    // outTree << maxLevel << endl;
    //
    // for (int i = 0; i < tree.size(); i++) {
    //     outTree << tree[i].name << endl;
    //     outTree << tree[i].level << '\t' << tree[i].parent << endl;
    // }
    //
    // //print genus nodes
    // outTree << endl << uniqueTaxonomies.size() << endl;
    // std::set<int>::iterator it2;
    // for (it2=uniqueTaxonomies.begin(); it2!=uniqueTaxonomies.end(); it2++) {  outTree << *it2 << '\t' << tree[*it2].accessions.size() << endl;	}
    // outTree << endl;
    //
    // outTree.close();
}
/**************************************************************************************************/
TaxNode PhyloTree::get(int i ){
	if (i < tree.size()) {  return tree[i];	 }
//	else {  m->mothurOut("Mismatch with taxonomy and template files. Cannot continue.\n");  exit(1); }
}
/**************************************************************************************************/
TaxNode PhyloTree::get(std::string seqName){
	std::map<std::string, int>::iterator itFind = name2Taxonomy.find(seqName);

	if (itFind != name2Taxonomy.end()) {  return tree[name2Taxonomy[seqName]];  }
	//else { m->mothurOut("Cannot find " + seqName + ". Mismatch with taxonomy and template files. Cannot continue.\n");  exit(1);}
}
/**************************************************************************************************/
std::string PhyloTree::getName(int i ){
	if (i < tree.size()) {  return tree[i].name;	 }
	//else { m->mothurOut("Mismatch with taxonomy and template files. Cannot continue.\n");  exit(1); }
}
/**************************************************************************************************/
int PhyloTree::getGenusIndex(std::string seqName){

	std::map<std::string, int>::iterator itFind = name2GenusNodeIndex.find(seqName);

	if (itFind != name2GenusNodeIndex.end()) {  return itFind->second;  }
	//else { m->mothurOut("Cannot find " + seqName + ". Could be a mismatch with taxonomy and template files. Cannot continue.\n");  exit(1);}
}
/**************************************************************************************************/
bool PhyloTree::ErrorCheck(std::vector<std::string> templateFileNames){
	bool okay = true;
	templateFileNames.push_back("unknown");

	std::map<std::string, int>::iterator itFind;
	std::map<std::string, int> taxonomyFileNames = name2Taxonomy;

  //  if (m->getDebug()) { m->mothurOut("[DEBUG]: in error check. Numseqs in template = " + tostd::string(templateFileNames.size()) + ". Numseqs in taxonomy = " + tostd::string(taxonomyFileNames.size()) + ".\n"); }

	for (int i = 0; i < templateFileNames.size(); i++) {
		itFind = taxonomyFileNames.find(templateFileNames[i]);

		if (itFind != taxonomyFileNames.end()) { //found it so erase it
			taxonomyFileNames.erase(itFind);
		}else {
		//	m->mothurOut("'" +templateFileNames[i] + "' is in your template file and is not in your taxonomy file. Please correct.\n");
			okay = false;
		}

		//templateFileNames.erase(templateFileNames.begin()+i);
		//i--;
	}
	templateFileNames.clear();

	if (taxonomyFileNames.size() > 0) { //there are names in tax file that are not in template
		okay = false;

		for (itFind = taxonomyFileNames.begin(); itFind != taxonomyFileNames.end(); itFind++) {
			//m->mothurOut(itFind->first + " is in your taxonomy file and is not in your template file. Please correct.\n");
		}
	}

	return okay;
}
/**************************************************************************************************/




