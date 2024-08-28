//
// Created by gregj on 8/27/2024.
//

#include "Adapters/CountTableAdapter.h"

void CountTableAdapter::CreateDataFrameMap(const Rcpp::DataFrame &countTable) {
    std::unordered_map<std::string, std::vector<double> > data;
    const std::vector<std::string> columnNames = countTable.names();
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
}

double CountTableAdapter::FindAbundanceBasedOnGroup(const std::string &group, const std::string &sampleName) const {
    if (std::find(groups.begin(), groups.end(), group) == groups.end())
        return -1; //Not Found, may need to throw and execption...
    const std::vector<double> groupCol = GetColumnByName(group);
    // Binary search?
    const long index = std::distance(sampleNames.begin(), std::find(sampleNames.begin(),
        sampleNames.end(), sampleName));
    return dataFrameMap.at(group)[index];
}

double CountTableAdapter::FindTotalAbundance(const std::string &sampleName) const{
    const long index = std::distance(sampleNames.begin(), std::find(sampleNames.begin(),
      sampleNames.end(), sampleName));
    return dataFrameMap.at("total")[index];
}

std::vector<std::string> CountTableAdapter::GetGroups() const {
    return groups;
}

std::vector<double> CountTableAdapter::GetColumnByName(const std::string &name) const {
    if (dataFrameMap.find(name) != dataFrameMap.end())
        return dataFrameMap.at(name);
    return {};
}
