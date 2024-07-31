//
// Created by Gregory Johnson on 7/31/24.
//

#include "SharedFileBuilder.h"
//TODO: Create Similar Sequences map **WE ARE HERE**
//TODO: Create a NameFile Based on the Similar Sequence Map
//TODO: Create a Counttable Based on the Name File
//TODO: Using the CountTable and a list file (obtained by clustering) create a Shared File
SharedFileBuilder::SharedFileBuilder(const RSpraseMatrix &sparseMatrix):
matrix(RSpraseMatrix(sparseMatrix)) {
}

void SharedFileBuilder::CreateNameMap() {
    // Creates the map of similar sequences neccessary for the name file!
    std::unordered_map<int, std::vector<int>> similarSequenes;
    for(const auto &row : matrix.sparseValues) {
        if(row.data < cutoff) {
            similarSequenes[row.iIndex].emplace_back(row.jIndex);
        }
    }
}
