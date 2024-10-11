//
// Created by Gregory Johnson on 10/3/24.
//

#include "Adapters/DistanceFilesToSparse.h"

#include <utility>

#include "MothurDependencies/ColumnDistanceMatrixReader.h"
#include "MothurDependencies/ReadPhylipMatrix.h"

DistanceFilesToSparse::DistanceFilesToSparse(CountTableAdapter countTable, std::string method, std::string filePath):
filePath(std::move(filePath)),
method(std::move(method)),
countTable(std::move(countTable)){}

std::vector<std::vector<double>> DistanceFilesToSparse::CreateSparseList() const {
    std::vector<RowData> data;
    if(method == "column") {
        ColumnDistanceMatrixReader reader(1.1, false);
        data = reader.ReadToRowData(countTable, filePath);
    }
    else {
        //phylip
        ReadPhylipMatrix reader(1.1, false); // We want every value
        data = reader.ReadToRowData(filePath);
    }
    const size_t size = data[0].rowValues.size();
    std::vector<double> iValues;
    std::vector<double> jValues;
    std::vector<double> dataValues;
    // order by index
    // Have to order them by there index relative to the count table
    const std::vector<std::string> samples = countTable.GetSamples();
    std::unordered_map<std::string, double> sampleToIndexMap;
    for(size_t i = 0; i < size; i++) {
        sampleToIndexMap[samples[i]] = static_cast<double>(i);
    }
   for(size_t i = 0; i < size; i++) {
       const std::string firstName = data[i].name;
       for(size_t j = 0; j < i; j++) {
           const double value = data[i].rowValues[j];
           if(value <= 0) continue;
           const std::string secondName = data[j].name;
           iValues.emplace_back(sampleToIndexMap[firstName]);
           jValues.emplace_back(sampleToIndexMap[secondName]);
           dataValues.emplace_back(value);
       }
   }
    return {iValues, jValues, dataValues};
}
