//
// Created by gregj on 8/27/2024.
//

#include "Adapters/CountTableAdapter.h"
#include "DataStructures/IndexAbundancePair.h"
#include "MothurDependencies/Utils.h"
#include <queue>

bool CountTableAdapter::CreateDataFrameMap(const Rcpp::DataFrame& count) {

    const std::vector<std::string> columnNames = count.names();

    if(columnNames.size() < 3)
        return false;
    bool sequenceColumn = true;
    for (const auto &name: columnNames) {
        if (sequenceColumn) {
            sequenceColumn = false; //Skip the first column
            // Sorted in R before we create the counttable adapter
            const std::vector<std::string> sequences = count[name];
            sequenceNames = sequences;
            continue;
        }
        // dataFrameMap[name] = Rcpp::as<std::vector<double>>(count[name]);
        dataFrameMap.emplace(name, Rcpp::as<std::vector<double>>(count[name]));
    }
    // In a count table, the first to columns are the sequence and the total abundance.
    // We only want the actual group names. so everything after
    groups.insert(groups.end(), columnNames.begin() + 2, columnNames.end());
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
            sequenceNames = samples;
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


        const size_t size = columnData.size();
        for(size_t i = 0; i < size; i++) {
            const std::string& columnString = columnData[i];
            if(columnString.empty()) {
                continue;
            }
            std::vector<std::string> container;
            Utils::splitAtComma(columnString, container);
            // the first one is the value we are looking for
            if(container[0] == "NA") continue;
            const int groupIndex = std::stoi(container[0]) - 1;
            const double abundance = std::stod(container[1]);
            // Since groupIndexes do not start until after the second column, we add two to the group index
            indexAbundanceQueue.push({groupIndex , static_cast<int>(i), abundance});
        }
    }
    groups.insert(groups.end(), columnNames.begin() + 2, columnNames.end());
    while(!indexAbundanceQueue.empty()) {
        const auto [groupIndex, sequenceIndex, abundance] = indexAbundanceQueue.front();
        data[groups[groupIndex]][sequenceIndex] = abundance;
        indexAbundanceQueue.pop();
    }
    dataFrameMap = data;
    // In a count table, the first to columns are the sequence and the total abundance.
    // We only want the actual group names. so everything after
    CreateNameToIndex();
    return true;

}

double CountTableAdapter::FindAbundanceBasedOnGroup(const std::string &group, const std::string &sampleName) const {
    // We will preprocess the find during hte read dist process. So remove special checks
    // - Protip hashmap find is faster than vector
    // if(nameToRowIndex.find(sampleName) == nameToRowIndex.end())
    //     return -1;
    return dataFrameMap.at(group)[nameToRowIndex.at(sampleName)];
}

double CountTableAdapter::FindTotalAbundance(const std::string &sampleName) const {
    // We will preprocess the find during hte read dist process. So remove special checks
    // - Protip hashmap find is faster than vector
    if(nameToRowIndex.find(sampleName) == nameToRowIndex.end())
        return -1;
    return dataFrameMap.at("total")[nameToRowIndex.at(sampleName)];
}

std::string CountTableAdapter::GetNameByIndex(const size_t index) const {
    if(index > sequenceNames.size())
        return "";
    return sequenceNames[index];
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
    for(const auto&[name, columnData]: dataFrameMap) {

        if(name == "total") {
            totals = columnData;
            continue;
        }
        const size_t index = groupIndexes[name];
        names[index] = name;
        columns[index] = columnData;
    }
    Rcpp::DataFrame countTable = Rcpp::DataFrame::create(Rcpp::Named("Representative Sequences") = sequenceNames,
        Rcpp::Named("total") = totals);
    for(size_t i = 0; i < size; i++) {
        countTable.push_back(columns[i], names[i]);
    }
    return countTable;
}

void CountTableAdapter::AddCountTable(const CountTableAdapter &other) {
    const std::vector<std::string>& sequences = other.GetSequences();
    const size_t currentSequenceSize = sequenceNames.size();
    sequenceNames.insert(sequenceNames.end(), sequences.begin(), sequences.end());
    const size_t otherSequenceSize = sequences.size();
    const size_t newSequenceSize = currentSequenceSize + otherSequenceSize;
    nameToRowIndex.reserve(newSequenceSize);

    for(size_t i = 0; i < otherSequenceSize; i++) {
        nameToRowIndex[sequences[i]] = i + currentSequenceSize;
    }

    groups.insert(groups.end(), other.groups.begin(), other.groups.end());
    for (const auto& group : groups) {
        const bool currentHasGroup = dataFrameMap.find(group) != dataFrameMap.end();
        const bool otherHasGroup = other.dataFrameMap.find(group) != other.dataFrameMap.end();

        if (otherHasGroup && currentHasGroup)  { // If they both have it
            const std::vector<double>& abundances = other.dataFrameMap.at(group);
            dataFrameMap[group].insert(dataFrameMap[group].end(), abundances.begin(),
                abundances.end());
            continue;
        }
        // It always means one has it and one doesnt have the group
        if (!otherHasGroup) {
            dataFrameMap[group].resize(newSequenceSize, 0);
            continue;
        }
        // Means !currentHasGroup
        std::vector<double> newAbundances = std::vector<double>(currentSequenceSize, 0);
        const std::vector<double>& abundances = other.dataFrameMap.at(group);
        newAbundances.insert(newAbundances.end(), abundances.begin(), abundances.end());
        dataFrameMap[group] = newAbundances;
    }

}


void CountTableAdapter::CreateNameToIndex() {
    for(size_t i = 0; i < sequenceNames.size(); i++) {
        nameToRowIndex[sequenceNames[i]] = i;
    }
}

// Gets every column but the first column (the sequence names)
std::vector<double> CountTableAdapter::GetColumnByName(const std::string &name) const {
    if (dataFrameMap.find(name) != dataFrameMap.end())
        return dataFrameMap.at(name);
    return {};
}
