//
// Created by Gregory Johnson on 9/2/26.
//

#include "ClusterSplit.h"

#include "../DataStructures/OptiData.h"
#include "../DataStructures/SplitMatrix.h"

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
    	// SplitMatrix* split = new SplitMatrix(taxaData, taxonomyCutoff);
        // m->mothurOut("Splitting the file...\n");
        // current->setMothurCalling(true);

        //split matrix into non-overlapping groups
        //SplitMatrix* split = new SplitMatrix(fastafile, namefile, countfile, taxFile, taxLevelCutoff, cutoff,  processors, classic, outputdir, usingVsearchToCLuster);

        //if (fastafile != "") {  current->setFastaFile(fastafile);  }

        // if (m->getControl_pressed()) { delete split; return 0; }

        // std::vector<std::string> singletonName = split->getSingletonNames();// create  multiple optimatrices
        // std::vector<std::map<std::string, std::string>> distName = split->getDistanceFiles();  //returns map of distance files -> namefile sorted by distance file size
        // delete split;
		// The split data
    	std::vector<OptiData*> splitMatrices =
    		SplitMatrix::splitClassify(taxaData, fastaData, calculator, cutoff, taxonomyCutoff);
        // current->setMothurCalling(false);

        // if (m->getDebug()) { m->mothurOut("[DEBUG]: distName.size() = " + std::to_string(distName.size()) + ".\n"); }

        Rcpp::message("It took " + std::to_string(time(nullptr) - estart) + " seconds to split the distance file.\n");

        // output a merged distance file
        // if (makeDist)		{ createMergedDistanceFile(distName); }

        // if (m->getControl_pressed()) { return 0; }

        estart = time(nullptr);
        //
        // if (!runCluster) {
        //     std::string filename = printFile(singletonName, distName);
        //
        //     m->mothurOutEndLine();
        //     m->mothurOut("Output File Names:\n\n"); m->mothurOut(filename); m->mothurOutEndLine();
        //     for (int i = 0; i < distName.size(); i++) {	m->mothurOut(distName[i].begin()->first); m->mothurOutEndLine(); m->mothurOut(distName[i].begin()->second); m->mothurOutEndLine();	}
        //     m->mothurOutEndLine();
        //
        //     return 0;
        // }
        bool deleteFiles = true;
    }
	//****************** break up files between processes and cluster each file set ******************************//

    listFileNames = createProcesses(distName, labels);

    // if (deleteFiles) {
    //     //delete the temp files now that we are done
    //     for (int i = 0; i < distName.size(); i++) {
    //         string thisNamefile = distName[i].begin()->second;
    //         string thisDistFile = distName[i].begin()->first;
    //         util.mothurRemove(thisNamefile);
    //         util.mothurRemove(thisDistFile);
    //     }
    // }

	// if (m->getControl_pressed()) { for (int i = 0; i < listFileNames.size(); i++) { util.mothurRemove(listFileNames[i]); } return 0; }

	if (!Utils::isEqual(saveCutoff, cutoff)) { Rcpp::message("\nCutoff was " + std::to_string(saveCutoff) +
		" changed cutoff to " + std::to_string(cutoff));  }

	Rcpp::message("It took " + std::to_string(time(nullptr) - estart) + " seconds to cluster\n");

	//****************** merge list file and create rabund and sabund files ******************************//
	estart = time(nullptr);
	//m->mothurOut("Merging the clustered files...\n");

	ListVector* listSingle;
	std::map<double, int> labelBins = completeListFile(listFileNames, singletonName, labels, listSingle); //returns map of label to numBins

//	if (m->getControl_pressed()) { if (listSingle != nullptr) { delete listSingle; } for (int i = 0; i < outputNames.size(); i++) { util.mothurRemove(outputNames[i]); } return 0; }

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
//**********************************************************************************************************************
std::map<double, int> ClusterSplit::completeListFile(std::vector<std::string> listNames, std::string singleton,
	std::set<std::string>& userLabels, ListVector*& listSingle){
		std::map<double, int> labelBin;
		std::vector<double> orderFloat;
		int numSingleBins;

		//read in singletons
		if (singleton != "none") {

            listSingle = new ListVector();
            if (type == "count") {

                CountTable ct; ct.readTable(singleton, false, false);

                std::vector<std::string> singletonSeqNames = ct.getNamesOfSeqs();

                for (int i = 0; i < singletonSeqNames.size(); i++) {
	                listSingle->push_back(singletonSeqNames[i]);
                }

            }else if (type == "name") {
                std::map<std::string, std::string> singletonSeqNames; Utils::readNames(singleton, singletonSeqNames);

                for (std::map<std::string, std::string>::iterator it = singletonSeqNames.begin(); it != singletonSeqNames.end(); it++) {
	                listSingle->push_back(it->second);
                }
            }

			util.mothurRemove(singleton);

			numSingleBins = listSingle->getNumBins();
        }else{  listSingle = nullptr; numSingleBins = 0;  }

        //go through users set and make them floats so we can sort them
        for(const auto & userLabel : userLabels) {
            double temp = -10.0;

            if ((userLabel != "unique") && (convertTestFloat(userLabel, temp) ))	{	util.mothurConvert(*it, temp);	}
            else if (userLabel == "unique")										{	temp = -1.0;		}

            if ((temp < cutoff) || util.isEqual(cutoff, temp)) {
                orderFloat.push_back(temp);
                labelBin[temp] = numSingleBins; //initialize numbins
            }
        }

		//sort order
		sort(orderFloat.begin(), orderFloat.end());
		userLabels.clear();

		//get the list info from each file
		for (int k = 0; k < listNames.size(); k++) {

			if (m->getControl_pressed()) {
				if (listSingle != nullptr) { delete listSingle; listSingle = nullptr; util.mothurRemove(singleton);  }
				for (int i = 0; i < listNames.size(); i++) {   util.mothurRemove(listNames[i]);  }
				return labelBin;
			}

			InputData* input = new InputData(listNames[k], "list", nullVector);
			ListVector* list = input->getListVector();
			string lastLabel = list->getLabel();

			string filledInList = listNames[k] + "filledInTemp";
			ofstream outFilled;
			util.openOutputFile(filledInList, outFilled);
            bool printHeaders = true;


			//for each label needed
			for(double & l : orderFloat){

				std::string thisLabel;
				if (Utils::isEqual(orderFloat[l],-1)) { thisLabel = "unique"; }
				else {
					thisLabel = std::to_string(l,  length-1);
				}

				//this file has reached the end
				if (list == nullptr) {
					list = input->getListVector(lastLabel, true);
				}else{	//do you have the distance, or do you need to fill in

					float labelFloat;
					if (const std::string& labelString = list->getLabel(); labelString == "unique") {
						labelFloat = -1.0;
					}
					else {
						labelFloat = atof(labelString);
						// convert(list->getLabel(), labelFloat);
					}

					//check for missing labels
					if (labelFloat > l) { //you are missing the label, get the next smallest one
						//if its bigger get last label, otherwise keep it
						delete list;
						list = input->getListVector(lastLabel, true);  //get last list vector to use, you actually want to move back in the file
					}
					lastLabel = list->getLabel();
				}

				//print to new file
				list->setLabel(thisLabel);
                list->setPrintedLabels(printHeaders);
                list->print(outFilled, true); printHeaders = false;

				//update labelBin
				labelBin[l] += list->getNumBins();

				delete list;

				list = input->getListVector();
			}

			if (list != nullptr) { delete list; }
			delete input;

			outFilled.close();
			util.mothurRemove(listNames[k]);
			rename(filledInList.c_str(), listNames[k].c_str());
		}

		return labelBin;
}


std::vector<std::string> ClusterSplit::createProcesses(std::vector<OptiData*>& distanceMatrices, std::set<std::string>& labels,
	size_t processors){
    //sanity check
	const size_t matricesSize = distanceMatrices.size();
    if (processors > matricesSize) {
	    processors = matricesSize;
    }
  //  deleteFiles = false; //so if we need to recalc the processors the files are still there
   // vector<string> listFiles;
    std::vector<std::vector<OptiData*>> dividedWork(processors); //distNames[1] = vector of filenames for process 1...
    // dividedNames.resize(processors);

    //for each file group figure out which process will complete it
    //want to divide the load intelligently so the big files are spread between processes
    for (int i = 0; i < matricesSize; i++) {
        int processToAssign = (i+1) % processors;
        if (processToAssign == 0) {
	        processToAssign = processors;
        }

        dividedWork[(processToAssign-1)].emplace_back(distanceMatrices[i]);
        // if ((processToAssign-1) == 1) { m->mothurOut(distName[i].begin()->first + "\n"); }
    }

    //now lets reverse the order of ever other process, so we balance big files running with little ones
    for (int i = 0; i < processors; i++) {
        int remainder = ((i+1) % processors);
        if (remainder) { std::reverse(dividedWork[i].begin(), dividedWork[i].end());  }
    }

    // if (m->getControl_pressed()) { return listFiles; }


    //create array of worker threads
    std::vector<std::thread*> workerThreads;
    std::vector<ClusterData*> data;

    //Lauch worker threads
    for (int i = 1; i < processors; i++) {
        // ClusterData* dataBundle = new ClusterData(showabund, classic, deleteFiles, dividedWork[i+1], cutoffNotSet, cutoff, precision, length, method, outputdir, vsearchLocation, type);
        // dataBundle->setOptiOptions(metricName, stableMetric, initialize, maxIters);

    	ClusterData* dataBundle = new ClusterData(dividedWork[i], );
        data.push_back(dataBundle);

        workerThreads.push_back(new std::thread(cluster, dataBundle));
    }


    ClusterData* dataBundle = new ClusterData(showabund, classic, deleteFiles, dividedWork[0], cutoffNotSet, cutoff, precision, length, method, outputdir, vsearchLocation, type);
    dataBundle->setOptiOptions(metricName, stableMetric, initialize, maxIters);
    cluster(dataBundle);
    listFiles = dataBundle->listFileNames;
    tag = dataBundle->tag;
    cutoff = dataBundle->cutoff;
    labels = dataBundle->labels;


    for (int i = 0; i < processors-1; i++) {
        workerThreads[i]->join();

        listFiles.insert(listFiles.end(), data[i]->listFileNames.begin(), data[i]->listFileNames.end());
        labels.insert(data[i]->labels.begin(), data[i]->labels.end());
        if (data[i]->cutoff < cutoff) { cutoff = data[i]->cutoff; }

        delete data[i];
        delete workerThreads[i];
    }
    delete dataBundle;
    deleteFiles = true;

    return listFiles;
}



void ClusterSplit::cluster(ClusterData* params){
	double smallestCutoff = params->cutoff;
	params->clusterMethod->Execute();
	//cluster each distance file
	for (int i = 0; i < params->distNames.size(); i++) {

		string thisNamefile = params->distNames[i].begin()->second;
		string thisDistFile = params->distNames[i].begin()->first;

		params->setNamesCount(thisNamefile);

		string listFileName = "";
		if (params->classic)    {  listFileName = clusterClassicFile(thisDistFile, thisNamefile, smallestCutoff, params);   }
		else                    {  listFileName = clusterFile(thisDistFile, thisNamefile, smallestCutoff, params);          }

		if (params->m->getControl_pressed()) { //clean up
			for (int i = 0; i < listFileNames.size(); i++) {	params->util.mothurRemove(listFileNames[i]); 	}
			params->listFileNames.clear(); break;
		}
		params->listFileNames.push_back(listFileName);
	}
	params->cutoff = smallestCutoff;
}

}