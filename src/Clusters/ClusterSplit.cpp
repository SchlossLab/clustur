//
// Created by Gregory Johnson on 9/2/26.
//

#include "ClusterSplit.h"

ClusterExport * ClusterSplit::Execute() {
    time_t estart;
    // vector<string> listFileNames;
    // vector< map<string, string> > distName;
    // set<string> labels;
    // string singletonName = "";

    double saveCutoff = 0;

    if (false) {
        // deleteFiles = false; estart = time(nullptr);
        // singletonName = readFile(distName);
        //
        // if (isList) {
        //
        //     //set list file as new current listfile
        //     string currentName = "";
        //     itTypes = outputTypes.find("list");
        //     if (itTypes != outputTypes.end()) {
        //         if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setListFile(currentName); }
        //     }
        //
        //     m->mothurOut("\nOutput File Names: \n");
        //     for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] + "\n");	} m->mothurOutEndLine();
        //
        //     return 0;
        // }

    }else {
        //splitting
        estart = time(nullptr); bool usingVsearchToCLuster = false;
        // if ((method == "agc") || (method == "dgc")) { usingVsearchToCLuster = true; if (cutoffNotSet) {  m->mothurOut("\nYou did not set a cutoff, using 0.03.\n"); cutoff = 0.03; } }

        // m->mothurOut("Splitting the file...\n");
        // current->setMothurCalling(true);

        //split matrix into non-overlapping groups
        SplitMatrix* split = new SplitMatrix(fastafile, namefile, countfile, taxFile, taxLevelCutoff, cutoff,  processors, classic, outputdir, usingVsearchToCLuster);

        if (fastafile != "") {  current->setFastaFile(fastafile);  }

        if (m->getControl_pressed()) { delete split; return 0; }

        singletonName = split->getSingletonNames();
        distName = split->getDistanceFiles();  //returns map of distance files -> namefile sorted by distance file size
        delete split;
        current->setMothurCalling(false);

        if (m->getDebug()) { m->mothurOut("[DEBUG]: distName.size() = " + toString(distName.size()) + ".\n"); }

        m->mothurOut("It took " + toString(time(nullptr) - estart) + " seconds to split the distance file.\n");

        //output a merged distance file
        if (makeDist)		{ createMergedDistanceFile(distName); }

        if (m->getControl_pressed()) { return 0; }

        estart = time(nullptr);

        if (!runCluster) {
            string filename = printFile(singletonName, distName);

            m->mothurOutEndLine();
            m->mothurOut("Output File Names:\n\n"); m->mothurOut(filename); m->mothurOutEndLine();
            for (int i = 0; i < distName.size(); i++) {	m->mothurOut(distName[i].begin()->first); m->mothurOutEndLine(); m->mothurOut(distName[i].begin()->second); m->mothurOutEndLine();	}
            m->mothurOutEndLine();

            return 0;
        }
        deleteFiles = true;
    }
	//****************** break up files between processes and cluster each file set ******************************//

    listFileNames = createProcesses(distName, labels);

    if (deleteFiles) {
        //delete the temp files now that we are done
        for (int i = 0; i < distName.size(); i++) {
            string thisNamefile = distName[i].begin()->second;
            string thisDistFile = distName[i].begin()->first;
            util.mothurRemove(thisNamefile);
            util.mothurRemove(thisDistFile);
        }
    }

	if (m->getControl_pressed()) { for (int i = 0; i < listFileNames.size(); i++) { util.mothurRemove(listFileNames[i]); } return 0; }

	if (!util.isEqual(saveCutoff, cutoff)) { m->mothurOut("\nCutoff was " + toString(saveCutoff) + " changed cutoff to " + toString(cutoff)); m->mothurOutEndLine();  }

	m->mothurOut("It took " + toString(time(nullptr) - estart) + " seconds to cluster\n");

	//****************** merge list file and create rabund and sabund files ******************************//
	estart = time(nullptr);
	m->mothurOut("Merging the clustered files...\n");

	ListVector* listSingle;
	map<double, int> labelBins = completeListFile(listFileNames, singletonName, labels, listSingle); //returns map of label to numBins

	if (m->getControl_pressed()) { if (listSingle != nullptr) { delete listSingle; } for (int i = 0; i < outputNames.size(); i++) { util.mothurRemove(outputNames[i]); } return 0; }

	mergeLists(listFileNames, labelBins, listSingle);

	if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) { util.mothurRemove(outputNames[i]); } return 0; }

    //delete after all are complete incase a crash happens
    if (!deleteFiles) { for (int i = 0; i < distName.size(); i++) {	util.mothurRemove(distName[i].begin()->first); util.mothurRemove(distName[i].begin()->second); 	} }

	m->mothurOut("It took " + toString(time(nullptr) - estart) + " seconds to merge.\n");

    if ((method == "opti") && (runsensSpec)) { runSensSpec();  }

    if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) { util.mothurRemove(outputNames[i]); } return 0; }

	//set list file as new current listfile
	string currentName = "";
	itTypes = outputTypes.find("list");
	if (itTypes != outputTypes.end()) {
		if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setListFile(currentName); }
	}

	//set rabund file as new current rabundfile
	itTypes = outputTypes.find("rabund");
	if (itTypes != outputTypes.end()) { if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setRabundFile(currentName); } }

	//set sabund file as new current sabundfile
	itTypes = outputTypes.find("sabund");
	if (itTypes != outputTypes.end()) { if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setSabundFile(currentName); } }

    //set sabund file as new current sabundfile
    itTypes = outputTypes.find("column");
    if (itTypes != outputTypes.end()) { if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setColumnFile(currentName); } }

	m->mothurOut("\nOutput File Names: \n");
	for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();

	return 0;
}
