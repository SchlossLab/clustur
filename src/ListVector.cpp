//
// Created by Gregory Johnson on 6/9/26.
//


//
// Created by Gregory Johnson on 3/29/24.
//
#include "DataStructures/ListVector.h"
#include "DataStructures/RAbundVector.h"
#include "Adapters/DataFrameAdapter.h"

std::string ListVector::getOTUName(const long long bin) {
    if (static_cast<long long>(binLabels.size()) > bin) {
    } else { getLabels(); }
    return binLabels[bin];
}

ListVector::ListVector(const ListVector& other) {
    for(const auto& currentData : other.data) {
        push_back(currentData);
    }
    label = other.label;
    printListHeaders = other.printListHeaders;
    binLabels = other.binLabels;
    otuTag = other.otuTag;
    numBins = other.numBins;
    numSeqs = other.numSeqs;
    maxRank = other.maxRank;
}

/***********************************************************************/

void ListVector::push_back(const std::string& seqNames) {
    data.push_back(seqNames);
    const int nNames = Utils::getNumNames(seqNames);

    numBins++;

    if (nNames > maxRank) { maxRank = nNames; }

    numSeqs += nNames;
}

void ListVector::push_back(const ListVector& listVector) {
    const int currentSize = size();
    const int newSize = currentSize + listVector.size();
    resize(newSize);
    for (int i = currentSize; i < newSize; i++) {
        set(i, listVector.get(i - currentSize));
    }
}

void ListVector::set(const int binNumber, const std::string &seqNames) {
    const int nNames_old = Utils::getNumNames(data[binNumber]);
    data[binNumber] = seqNames;
    const int nNames_new = Utils::getNumNames(seqNames);

    if(nNames_old == 0)			{	numBins++;				}
    if(nNames_new == 0)			{	numBins--;				}
    if(nNames_new > maxRank)	{	maxRank = nNames_new;	}

    numSeqs += (nNames_new - nNames_old);
}

std::vector<std::string> ListVector::getLabels() {
    Utils::getOTUNames(binLabels, numBins, otuTag);
    return binLabels;
}

std::string ListVector::get(const long long index) const {
    if (index >= 0 && index < static_cast<long long>(data.size())) { return data[index]; }
    return "";
}

void ListVector::setLabels(const std::vector<std::string> &labels) {
    binLabels = labels;
    getLabels();
}

std::string ListVector::print(std::ostream &output, std::map<std::string, int> &ct) {
    otuTag = "Otu";
    std::string output_cluster;
    printHeaders(output_cluster, ct, true);
    std::vector<listCt> hold;
    for (const auto & i : data) {
        if (!i.empty()) {
            std::vector<std::string> binNames;
            std::string bin = i;
            Utils::splitAtComma(bin, binNames);
            int total = 0;
            for (const auto & binName : binNames) {
                auto it = ct.find(binName);
                if (it == ct.end()) {
                } else { total += ct[binName]; }
            }
            listCt temp(i, total, "");
            hold.push_back(temp);
        }
    }
    std::sort(hold.begin(), hold.end(), abundNamesSort2);
    for (auto & i : hold) {
        if (!i.bin.empty()) {
            //TestHelper::Print('\t' + hold[i].bin);
            output_cluster += "\t" + i.bin;
        }
    }
    output_cluster += "\n";
    return output_cluster;
}

int ListVector::size() const {
    return static_cast<int>(data.size());
}

void ListVector::clear() {
    data.clear();
}
void ListVector::resize(const int size) {
    data.resize(size);
}

std::string ListVector::print(std::ostream &output) {
    std::map<std::string, int> ct;
    for (const auto & i : data) {
        if (!i.empty()) {
            std::string bin = i;
            std::vector<std::string> binNames;
            Utils::splitAtComma(bin, binNames);
            if (!std::isdigit(bin[0])) {
                binNames[0] = "";
                //continue;
            }
            for (const auto& key : binNames) {
                ct[key] = 1;
            }
        }
    }

    return print(output, ct);
}
RAbundVector ListVector::getRAbundVector() const {
    RAbundVector rav;
    for(const auto & i : data){
        const int binSize = Utils::getNumNames(i);
        rav.push_back(binSize);
    }
    rav.setLabel(label);

    return rav;
}

void ListVector::printHeaders(std::string &output, std::map<std::string, int> &ct, bool sortPlease) {
    if (printListHeaders) {
        if (binLabels.empty()) { sortPlease = false; } //we are creating arbitary otuNames
        const std::vector<std::string> theseLabels = getLabels();
        if (sortPlease) {
            std::vector<listCt> hold;
            for (size_t i = 0; i < data.size(); i++) {
                if (!data[i].empty()) {
                    std::vector<std::string> binNames;
                    std::string bin = data[i];
                    Utils::splitAtComma(bin, binNames);
                    int total = 0;
                    for (const auto & binName : binNames) {
                        auto it = ct.find(binName);
                        if (it == ct.end()) {
                            //m->mothurOut("[ERROR]: " + binNames[j] + " is not in your count table. Please correct.\n"); m->setControl_pressed(true);
                        } else { total += ct[binName]; }
                    }
                    listCt temp(data[i], total, theseLabels[i]);
                    hold.push_back(temp);
                }
            }
            std::sort(hold.begin(), hold.end(), abundNamesSort2);

            //print original label for sorted by abundance otu
            for (auto & i : hold) { output += ('\t' + i.label); }
        } else {
            for (const auto & theseLabel : theseLabels) { output += ('\t' + theseLabel); }
        }

        output += "\n";
        printListHeaders = false;
    }
}
///Create Test
Rcpp::DataFrame ListVector::CreateDataFrameFromList(const std::string& featureColumnName,
    const std::string& binColumnName) const {
    std::unordered_map<std::string, std::vector<std::string>> map;
    const std::vector<std::string> headers{binColumnName, featureColumnName};
    int count = 1;
    for(const auto& bin : data) {
        if(bin.empty())
            continue;
        map[headers[0]].emplace_back(binColumnName + std::to_string(count++));
        map[headers[1]].emplace_back(bin);
    }
    return DataFrameAdapter::UnorderedMapToDataFrame(map);
}
