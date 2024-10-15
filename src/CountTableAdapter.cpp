//
// Created by gregj on 8/27/2024.
//

#include "Adapters/CountTableAdapter.h"

#include "MothurDependencies/Utils.h"

bool CountTableAdapter::CreateDataFrameMap(const Rcpp::DataFrame &countTable) {
    
    std::unordered_map<std::string, std::vector<double> > data;
    const std::vector<std::string> columnNames = countTable.names();

    if(columnNames.size() < 3)
        return false;
    bool sequenceColumn = true;
    for (const auto &name: columnNames) {
        if (sequenceColumn) {
            sequenceColumn = false; //Skip the first column
            const std::vector<std::string> samples = countTable[name];
            sampleNames = samples;
            continue;
        }
        const std::vector<double> columnData = countTable[name];
        data[name] = columnData;
    }
    dataFrameMap = data;
    // In a count table, the first to columns are the sequence and the total abundance.
    // We only want the actual group names. so everything after
    groups.insert(groups.end(), columnNames.begin() + 2, columnNames.end());
    this->countTable = countTable;
    CreateNameToIndex();
    return true;
}

bool CountTableAdapter::CreateDataFrameMapFromSparseCountTable(const Rcpp::DataFrame &countTable) {
    std::unordered_map<std::string, std::vector<double> > data;
    std::queue<IndexAbundancePair> indexAbundanceQueue;
    const std::vector<std::string> columnNames = countTable.names();
    const auto rowSize = static_cast<size_t>(countTable.nrows());
    if(columnNames.size() < 3)
        return false;

    bool sequenceColumn = true;
    for (const auto &name: columnNames) {
        if (sequenceColumn) {
            sequenceColumn = false; //Skip the first column
            const std::vector<std::string> samples = countTable[name];
            sampleNames = samples;
            continue;
        }
        data[name] = std::vector<double>(rowSize, 0);
        const std::vector<std::string> columnData = countTable[name];
        if(name == "total") {
            const size_t size = columnData.size();
            std::vector<double> totals(size);
            for(size_t i = 0; i <  size; i++) {
                totals[i] = std::stod(columnData[i]);
            }
            data[name] = totals;
            continue;
        }


        Utils util;
        const size_t size = columnData.size();
        for(size_t i = 0; i < size; i++) {
            const std::string& columnString = columnData[i];
            if(columnString.empty()) {
                continue;
            }
            std::vector<std::string> container;
            util.splitAtComma(columnString, container);
            // the first one is the value we are looking for
            const int groupIndex = std::stoi(container[0]) - 1;
            const double abundance = std::stod(container[1]);
            // Since groupIndexes do not start until after the second column, we add two to the group index
            indexAbundanceQueue.push({groupIndex , static_cast<int>(i), abundance});
        }
    }
    groups.insert(groups.end(), columnNames.begin() + 2, columnNames.end());
    while(!indexAbundanceQueue.empty()) {
        const IndexAbundancePair pair = indexAbundanceQueue.front();
        data[groups[pair.groupIndex]][pair.sequenceIndex] = pair.abundance;
        indexAbundanceQueue.pop();
    }
    dataFrameMap = data;
    // In a count table, the first to columns are the sequence and the total abundance.
    // We only want the actual group names. so everything after
    this->countTable = countTable;
    CreateNameToIndex();
    return true;

}

double CountTableAdapter::FindAbundanceBasedOnGroup(const std::string &group, const std::string &sampleName) const {
    // We will preprocess the find during hte read dist process. So remove special checks
    // - Protip hashmap find is faster than vector
    if(nameToRowIndex.find(sampleName) == nameToRowIndex.end())
        return -1;
    const std::vector<double> groupCol = GetColumnByName(group);
    return dataFrameMap.at(group)[nameToRowIndex.at(sampleName)];
}

double CountTableAdapter::FindTotalAbundance(const std::string &sampleName) const {
    // We will preprocess the find during hte read dist process. So remove special checks
    // - Protip hashmap find is faster than vector
    if(nameToRowIndex.find(sampleName) == nameToRowIndex.end())
        return -1;
    return dataFrameMap.at("total")[nameToRowIndex.at(sampleName)];
}

std::string CountTableAdapter::GetNameByIndex(const int index) const {
    if(index > static_cast<int>(sampleNames.size()))
        return "";
    return sampleNames[index];
}

std::vector<std::string> CountTableAdapter::GetGroups() const {
    return groups;
}

Rcpp::DataFrame CountTableAdapter::ReCreateDataFrame() const {
    const size_t size = dataFrameMap.size() - 1;
    std::vector<double> totals;
    std::vector<std::vector<double>> columns(size);
    std::vector<std::string> names(size);
    std::unordered_map<std::string, size_t> groupIndexes;
    for(size_t i = 0; i < groups.size(); i++) {
        groupIndexes[groups[i]] = i;
    }
    for(const auto& column: dataFrameMap) {

        if(column.first == "total") {
            totals = column.second;
            continue;
        }
        const size_t index = groupIndexes[column.first];
        names[index] = column.first;
        columns[index] = column.second;
    }
    Rcpp::DataFrame countTable = Rcpp::DataFrame::create(Rcpp::Named("Representative Sequences") = sampleNames,
        Rcpp::Named("total") = totals);
    for(size_t i = 0; i < size; i++) {
        countTable.push_back(columns[i], names[i]);
    }
    return countTable;
}


void CountTableAdapter::CreateNameToIndex() {
    for(size_t i = 0; i < sampleNames.size(); i++) {
        nameToRowIndex[sampleNames[i]] = i;
    }
}

// Gets every column but the first column (the sequence names)
std::vector<double> CountTableAdapter::GetColumnByName(const std::string &name) const {
    if (dataFrameMap.find(name) != dataFrameMap.end())
        return dataFrameMap.at(name);
    return {};
}
