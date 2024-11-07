//
// Created by Gregory Johnson on 9/13/24.
//

#include "MothurDependencies/ColumnDistanceMatrixReader.h"
#include <Rcpp.h>

ColumnDistanceMatrixReader::ColumnDistanceMatrixReader(const double cutoff, const bool isSimularity)
:DistanceFileReader() {
	this->cutoff = cutoff;
	this->sim = isSimularity;
	sparseMatrix = new SparseDistanceMatrix();
	list = new ListVector();
}
bool ColumnDistanceMatrixReader::Read(const std::string& filePath) {
	std::ifstream fileHandle;
	fileHandle.open(filePath);
	if(!fileHandle.is_open())
		return false;
	Utils util;

	std::string firstName, secondName;
	float distance;
	std::vector<std::string> sequences = countTable.GetSamples();
	size_t nseqs = sequences.size();
    sparseMatrix->resize(nseqs);
	list = new ListVector(static_cast<int>(nseqs));
	std::unordered_map<std::string, int> nameToIndexMap;
	int count = 0;
	for(const auto &sequence : sequences) {
		list->set(count, sequence);
		nameToIndexMap[sequence] = count++;
	}

    int lt = 1;
	int refRow = 0;	//we'll keep track of one cell - Cell(refRow,refCol) - and see if it's transpose
	int refCol = 0; //shows up later - Cell(refCol,refRow).  If it does, then its a square matrix

	//need to see if this is a square or a triangular matrix...
	std::string dist;
	fileHandle >> firstName;
	fileHandle >> secondName;
	fileHandle >> dist;

	if(nameToIndexMap.find(firstName) != nameToIndexMap.end() ||
			nameToIndexMap.find(secondName) != nameToIndexMap.end()) {
			fileHandle.clear();
			fileHandle.seekg(0, std::ifstream::beg);
	}
	while(fileHandle >> firstName >> secondName >> distance && lt == 1){  //let's assume it's a triangular matrix...
        int itA;
		int itB;
		try {
			itA = nameToIndexMap.at(firstName);
			itB = nameToIndexMap.at(secondName);
		}
		catch (const std::exception& ex) {
			std::set<std::string> container;
			if(nameToIndexMap.find(firstName) == nameToIndexMap.end()) {
				container.insert(firstName);
			}
			if(nameToIndexMap.find(secondName) == nameToIndexMap.end()) {
				container.insert(secondName);
			}
			util.CheckForDistanceFileError(container);
		}

		if (util.isEqual(distance, -1)) { distance = 1000000; }
		else if (sim) { distance = 1 - distance;  }  //user has entered a sim matrix that we need to convert.

		if(distance <= cutoff && itA != itB){
			if(itA > itB){
                PDistCell value(itA, distance);


				if(refRow == refCol){		// in other words, if we haven't loaded refRow and refCol...
					refRow = itA;
					refCol = itB;
					sparseMatrix->addCell(itB, value); // This is the problem most likely...How do we fix it.
				}
				else if(refRow == itA && refCol == itB){
					lt = 0;
				}
				else{
					sparseMatrix->addCell(itB, value);
				}
			}
			else if(itA < itB){
				PDistCell value(itB, distance);

				if(refRow == refCol){		// in other words, if we haven't loaded refRow and refCol...
					refRow = itA;
					refCol = itB;
					sparseMatrix->addCell(itA, value);
				}
				else if(refRow == itB && refCol == itA){
					lt = 0;
				}
				else{
					sparseMatrix->addCell(itA, value);
				}
			}
		}
	}

	if(lt == 0) {  // oops, it was square
		fileHandle.close();  //let's start over
		sparseMatrix->clear();  //let's start over
		fileHandle.open(filePath); //let's start over

		while(fileHandle >> firstName >> secondName >> distance){
			int itA;
			int itB;
			try {
				itA = nameToIndexMap.at(firstName);
				itB = nameToIndexMap.at(secondName);
			}
			catch (const std::exception& ex) {
				std::set<std::string> container;
				if(nameToIndexMap.find(firstName) == nameToIndexMap.end()) {
					container.insert(firstName);
				}
				if(nameToIndexMap.find(secondName) == nameToIndexMap.end()) {
					container.insert(secondName);
				}
				util.CheckForDistanceFileError(container);
			}

			if (util.isEqual(distance, -1)) { distance = 1000000; }
			else if (sim) { distance = 1 - distance;  }  //user has entered a sim matrix that we need to convert.

			if(distance <= cutoff && itA > itB){
                PDistCell value(itA, distance);
				sparseMatrix->addCell(itB, value);
			}
		}
	}
	fileHandle.close();
	list->setLabel("0");
	return true;
}


