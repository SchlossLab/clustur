//
// Created by Gregory Johnson on 7/31/24.
//

#include "MothurDependencies/SharedFileBuilder.h"


//TODO: Create Similar Sequences map **WE ARE HERE**
//TODO: Create a NameFile Based on the Similar Sequence Map
//TODO: Create a Counttable Based on the Name File
//TODO: Using the CountTable and a list file (obtained by clustering) create a Shared File
SharedFileBuilder::SharedFileBuilder(const RSpraseMatrix &sparseMatrix):
matrix(RSpraseMatrix(sparseMatrix)) {
}

void SharedFileBuilder::CreateNameMap() {
    // Creates the map of similar sequences neccessary for the name file!
    std::unordered_map<int, std::set<int>> similarSequenes;
    for(const auto& row : matrix.sparseValues) {
        similarSequenes[row.iIndex].insert(row.iIndex);
    }

    for(const auto &row : matrix.sparseValues) {
        if(row.data < cutoff) {
            similarSequenes[row.iIndex].insert(row.jIndex);
        }
    }
    nameMap = similarSequenes;
}

void SharedFileBuilder::CreateCountTable() {
    std::unordered_map<int, int> counts;
    for(const auto& row : matrix.sparseValues) {
        counts[row.iIndex] = static_cast<int>(nameMap[row.iIndex].size());
    }
    countTable = counts;

}

void SharedFileBuilder::OutputFiles() const {
    std::string nameFileData =  "Representative_Sequence\ttotal\n";
    std::string countTableFileData = "Sequence\tSimilar_Sequences\n";
    for(const auto& values : nameMap) {
        nameFileData += (std::to_string(values.first) + "\t");
        std::string otherNames;
        for(auto const& names : values.second) {
            nameFileData += (std::to_string(names) + ",");
        }
        nameFileData[nameFileData.size() - 1] = '\n';
    }
    for(const auto& values : countTable) {
        countTableFileData += std::to_string(values.first) + "\t" +
            std::to_string(values.second) + "\n";
    }
    std::ofstream nameFileOutput("/Users/grejoh/Documents/OptiClusterPackage/Opticluster/NameFile.txt");
    if(!nameFileOutput.is_open()) {
        std::cout << "Please enter a new name file name";
    }
    nameFileOutput.close();
    std::ofstream countTableOutput("/Users/grejoh/Documents/OptiClusterPackage/Opticluster/CountTable.txt");
    nameFileOutput << nameFileData;
    countTableOutput << countTableFileData;
    nameFileOutput.close();
}

SharedFile* SharedFileBuilder::BuildSharedFile(ListVector *listVector) {
    const int size = listVector->size();
    std::vector<OTUAbundancePair> abundanceMap(size);
    std::vector<OTUAbundancePair> pair;
    for(int i = 0; i < size; i++) {
        std::vector<std::string> container;
        std::string bin = listVector->get(i);
        utils.splitAtComma(bin, container);
        int abundance = 0;
        for(const auto& sequence : container) {
            abundance += countTable[std::stoi(sequence)];
        }
        abundanceMap[i] = OTUAbundancePair(abundance,
            "OTU" + std::to_string(i + 1),"");
        // Get abundance at each otu level
    }
    return new SharedFile(abundanceMap);
}
