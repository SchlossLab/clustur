//
// Created by Gregory Johnson on 6/9/26.
//

#ifndef REFACTOR_OPTIDATA_H
#define REFACTOR_OPTIDATA_H
#include <vector>
#include <string>
#include <unordered_set>
#include "../DataStructures/ListVector.h"
#include "../MothurDependencies/PairwiseDistanceCalculator.h"


class OptiData {
public:
    OptiData() = default;
    OptiData(const FastaDatabase& database, const PairwiseDistanceCalculator* calculator,
        const std::unordered_set<std::string>& filteredNameList, const double cutoff):cutoff(cutoff) {
        std::unordered_map<std::string, size_t> nameToIndex;
        nameToIndex.reserve(filteredNameList.size());
        const std::vector<FastaData>& data = database.GetFastaDataBase();
        nameMap.resize(filteredNameList.size());
        size_t count = 0;
        for (const auto& name : filteredNameList) {
            nameMap[count] = name;
            nameToIndex[name] = count++;
        }
        for (const auto& [name, sequence] : data) {
            if (filteredNameList.find(name) == filteredNameList.end()) continue;
            const size_t& index = nameToIndex[name];
            for (const auto& [otherName, otherSequence] : data) {
                if (filteredNameList.find(name) == filteredNameList.end()) continue;
                if (calculator->Execute(name, otherName) > cutoff) continue;
                size_t otherIndex = nameToIndex[otherName];
                closeness[index].emplace_back(otherIndex);
                closeness[otherIndex].emplace_back(index);
            }
        }
        for (const auto& [name, index] : nameToIndex) {
            if (closeness[index].size() < 0) {
                singletons.emplace_back(name);
            }
        }
    }
    virtual ~OptiData() = default;
    void MoveData(OptiData* other) {
        cutoff = other->cutoff;
        closeness = std::move(other->closeness);
        singletons = std::move(other->singletons);
        nameMap = std::move(other->nameMap);
    }
    [[nodiscard]] std::vector<std::vector<long long>> GetCloseness() const {return closeness;}
    [[nodiscard]] std::vector<std::string> GetNameList() const {return nameMap;}
    [[nodiscard]] std::vector<std::string> GetSingletons() const {return singletons;}
    [[nodiscard]] size_t getNumSeqs() const { return closeness.size(); }
    [[nodiscard]] size_t getNumSingletons() const { return singletons.size(); }
    [[nodiscard]] double GetCutoff() const {return cutoff;}
    [[nodiscard]] virtual bool isCloseFit(long long j, long long i, bool& k) const { return false; }
    virtual std::vector<long long> getTranslatedBins(std::vector<std::vector<std::string> >&,
        std::vector< std::vector<long long> >&) { return {};  }

    [[nodiscard]] bool isClose(const long long i, const long long toFind) const {
        if (i < 0) { return false; }
        if (i > static_cast<long long>(closeness.size())) {
            return false;
        }
        return std::binary_search(closeness[i].begin(), closeness[i].end(), toFind);
        // bool found = false;
        // if (closeness[i].find(toFind) != closeness[i].end()) { found = true; }
        // return found;

    }
    [[nodiscard]] std::vector<long long> getCloseSeqs(const long long i) const {
        if (i < 0) {
            return {};
        }
        if (i > static_cast<long long>(closeness.size())) {
            return {};
        }
        return closeness[i];
    }

    [[nodiscard]] long long getNumClose(const long long index) const {
        if (index < 0)
            return 0;
        if (index > static_cast<long long>(closeness.size())) {
            return 0;
        }
        return closeness[index].size();

    }

    [[nodiscard]] ListVector* getListSingle() const {
        ListVector* singlelist = nullptr;
        if (singletons.empty()) { }
        else {
            singlelist = new ListVector();

            for (const auto & singleton : singletons) {
                singlelist->push_back(singleton);
            }
        }
        return singlelist;

    }

    [[nodiscard]] std::string getName(const size_t index) const {

        if (index > closeness.size()) {
            return "";
        }

        return nameMap.at(index);

    }
    [[nodiscard]] bool mccValidCalc() const {

        const auto numSeqs = static_cast<double>(getNumSeqs());
        const double numDists = numSeqs * (numSeqs-1)/2;

        size_t totalClose = 0;
        //for each sequence (singletons removed on read)
        for (const auto & closenessData : closeness) {
            totalClose += closenessData.size();
        }

        // Inital setup of all singletons - badState <- TN == 0, FP == 0, FN == totalClose/2, TP = 0
        // Inital setup of one otu - badState <- TN == 0, FP == 0, FN == 0, TP = totalClose/2
        if ((numDists - (static_cast<double>(totalClose)/2)) == 0) {
            return false;
        }

        return true;
    }

    [[nodiscard]] virtual std::map<std::string, long long> getNameIndexMap() const {
        std::map<std::string, long long> nameIndexes;
        for (int i = 0; i < nameMap.size(); i++) {
            std::vector<std::string> thisBinsSeqs; Utils::splitAtComma(nameMap[i], thisBinsSeqs);
            if (i < closeness.size()) { nameIndexes[thisBinsSeqs[0]] = i;  }
        }
        return nameIndexes;
    }

    [[nodiscard]] virtual long long getNumDists() const {
        size_t foundDists = 0;
        for (const auto & closenes : closeness) {
            foundDists += closenes.size();
        }
        return static_cast<long long>(foundDists);
    }

protected:
    double cutoff{};
    std::vector<std::vector<long long>> closeness;  //closeness[0] contains indexes of seqs "close" to seq 0.
    std::vector<std::string> singletons; //name of seqs with NO distances in matrix, if name file is given then it contains 2nd column of namefile
    std::vector<std::string> nameMap;
};


#endif //REFACTOR_OPTIDATA_H
