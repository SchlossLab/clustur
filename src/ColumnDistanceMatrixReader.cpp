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
bool ColumnDistanceMatrixReader::Read(const CountTableAdapter &countTable, const std::string& filePath) {
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

	while(fileHandle && lt == 1){  //let's assume it's a triangular matrix...

		fileHandle >> firstName;
        fileHandle >> secondName;
        fileHandle >> distance;	// get the row and column names and distance
        int itA = nameToIndexMap[firstName];
		int itB = nameToIndexMap[secondName];

        // std::map<std::string,int>::iterator itB = nameMap->find(secondName);

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

	if(lt == 0){  // oops, it was square
		fileHandle.close();  //let's start over
		sparseMatrix->clear();  //let's start over
		fileHandle.open(filePath); //let's start over

		while(fileHandle){
			fileHandle >> firstName;
            fileHandle >> secondName;
            fileHandle >> distance;	// get the row and column names and distance


			int itA = nameToIndexMap[firstName];
			int itB = nameToIndexMap[secondName];

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


std::vector<RowData> ColumnDistanceMatrixReader::ReadToRowData(const CountTableAdapter &countTable,
	const std::string& filePath) {
	std::ifstream fileHandle;
	fileHandle.open(filePath);
	if(!fileHandle.is_open())
		return {};
	Utils util;

	std::string firstName, secondName;
	float distance;
	std::vector<std::string> sequences = countTable.GetSamples();
	size_t nseqs = sequences.size();
	std::vector<RowData> rowData(nseqs);
    sparseMatrix->resize(nseqs);
	list = new ListVector(static_cast<int>(nseqs));
	std::unordered_map<std::string, int> nameToIndexMap;
	std::map<std::string, int> name;
	for(const auto &sequence : sequences) {
		name[sequence] = 1;
	}
	int count = 0;
	for(const auto &sequence : name) {
		list->set(count, sequence.first);
		rowData[count].rowValues = std::vector<double>(nseqs, -1); // Values that surpass the threshold will have all 1s
		rowData[count].name = sequence.first;
		nameToIndexMap[sequence.first] = count++;
	}
    int lt = 1;
	int refRow = 0;	//we'll keep track of one cell - Cell(refRow,refCol) - and see if it's transpose
	int refCol = 0; //shows up later - Cell(refCol,refRow).  If it does, then its a square matrix

	//need to see if this is a square or a triangular matrix...

	while(fileHandle && lt == 1){  //let's assume it's a triangular matrix...

		fileHandle >> firstName;
        fileHandle >> secondName;
        fileHandle >> distance;	// get the row and column names and distance
        int itA = nameToIndexMap[firstName];
		int itB = nameToIndexMap[secondName];
        // std::map<std::string,int>::iterator itB = nameMap->find(secondName);

		if (util.isEqual(distance, -1)) { distance = 1000000; }
		else if (sim) { 
			distance = 1 - distance; 
		}  //user has entered a sim matrix that we need to convert.
		// if(distance > 0)
		// 	Rcpp::Rcout << "Distance: " << distance << std::endl;
		if(itA != itB){
			if(itA > itB){
                PDistCell value(itA, distance);
				if(refRow == refCol){		// in other words, if we haven't loaded refRow and refCol...
					refRow = itA;
					refCol = itB;
					rowData[itB].rowValues[itA] = distance;
					rowData[itA].rowValues[itB] = distance;
				}
				else if(refRow == itA && refCol == itB){
					lt = 0;
				}
				else{
					rowData[itB].rowValues[itA] = distance;
					rowData[itA].rowValues[itB] = distance;
				}
			}
			else if(itA < itB){
				if(refRow == refCol){		// in other words, if we haven't loaded refRow and refCol...
					refRow = itA;
					refCol = itB;
					rowData[itB].rowValues[itA] = distance;
					rowData[itA].rowValues[itB] = distance;
				}
				else if(refRow == itB && refCol == itA){
					lt = 0;
				}
				else{
					rowData[itB].rowValues[itA] = distance;
					rowData[itA].rowValues[itB] = distance;
				}
			}
		}
		else{rowData[itA].rowValues[itB] = 0;} // If itselfs, it should be zero
	}

	if(lt == 0){  // oops, it was square
		fileHandle.close();  //let's start over
		sparseMatrix->clear();  //let's start over
		fileHandle.open(filePath);
		firstName.clear();//let's start over
	// Clear rowData
		while(fileHandle){
			fileHandle >> firstName;
            fileHandle >> secondName;
            fileHandle >> distance;	// get the row and column names and distance

			int itA = nameToIndexMap[firstName];
			int itB = nameToIndexMap[secondName];

			if (util.isEqual(distance, -1)) {
				 distance = 1000000; 
			}
			else if (sim) { distance = 1 - distance;  }
			//if(distance > 0)
			//	Rcpp::Rcout << "Distance: " << distance << std::endl; //user has entered a sim matrix that we need to convert.
			rowData[itB].rowValues[itA] = distance;
			rowData[itA].rowValues[itB] = distance;
		}
	}
	fileHandle.close();
	list->setLabel("0");

	return rowData;

}



