//
// Created by Gregory Johnson on 9/13/24.
//

#include "ColumnDistanceMatrixReader.h"

#include <utility>

ColumnDistanceMatrixReader::ColumnDistanceMatrixReader(std::string path,
    const double cutoff):filePath(std::move(path)), cutoff(cutoff) {
}

bool ColumnDistanceMatrixReader::Read(const CountTableAdapter &countTable) {
	std::string firstName, secondName;
    std::ifstream fileHandle;
	Utils util;
    fileHandle.open(filePath);
    if(!fileHandle.is_open())
        return false;
	float distance;
	const std::vector<std::string> sequences = countTable.GetSamples();
    const size_t nseqs = sequences.size();
	std::unordered_map<std::string, int> nameToIndexMap;
    std::vector<RowData> rowData(nseqs);
	for(auto& row : rowData) {
		row.rowValues = std::vector<double>(nseqs);
	}
    sparseMatrix->resize(nseqs);
    list = new ListVector();
	int count = 0;
	for(const auto &sequence : sequences) {
		list->set(count, sequence);
		nameToIndexMap[sequence] = count++;

	}
	while(fileHandle){  //let's assume it's a triangular matrix...
		fileHandle >> firstName;
        fileHandle >> secondName;
        fileHandle >> distance;	// get the row and column names and distance

        int indexA = nameToIndexMap[firstName];
        int indexB = nameToIndexMap[secondName];
		if (util.isEqual(distance, -1)) { distance = 1000000; }
		else if (sim) { distance = 1 - distance;  }  //user has entered a sim matrix that we need to convert.

		if(distance <= cutoff && firstName != secondName){
			PDistCell value(indexA, distance);
			sparseMatrix->addCell(indexB, value);
		}
	}
	fileHandle.close();

	list->setLabel("0");

	return true;
}


std::vector<RowData> ColumnDistanceMatrixReader::readToRowData(const CountTableAdapter &countTable) {
    const size_t nseqs = countTable.GetColumnByName("total").size();
    std::vector<RowData> rowData(nseqs);
	return {};

}



