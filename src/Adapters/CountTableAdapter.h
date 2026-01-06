//
// Created by gregj on 8/27/2024.
//

#ifndef COUNTTABLEADAPTER_H
#define COUNTTABLEADAPTER_H
#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <queue>
#include <set>


class CountTableAdapter {
public:
    CountTableAdapter() = default;
    bool CreateDataFrameMap(const Rcpp::DataFrame& count);
    bool CreateDataFrameMapFromSparseCountTable(const Rcpp::DataFrame& countTable);
    std::vector<std::string> GetSequences() const {return sequenceNames;}
    // Going to ensure that each count_table atleast has a group
    // And if there is no count table inputted, I will create a base one.
    double FindAbundanceBasedOnGroup(const std::string& group, const std::string& sampleName) const;
    double FindTotalAbundance(const std::string& sampleName) const;
    std::string GetNameByIndex(size_t) const;
    std::vector<double> GetColumnByName (const std::string& name) const;
    std::vector<std::string> GetGroups() const;
    Rcpp::DataFrame ReCreateDataFrame() const;
private:
    void CreateNameToIndex();
    std::unordered_map<std::string, size_t> nameToRowIndex;
    std::vector<std::string> sequenceNames;
    std::unordered_map<std::string, Rcpp::NumericVector> dataFrameMap;
    std::vector<std::string> groups;
};



#endif //COUNTTABLEADAPTER_H
