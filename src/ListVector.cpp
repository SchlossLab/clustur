//
// Created by Gregory Johnson on 3/29/24.
//
#include "MothurDependencies/ListVector.h"
#include "Adapters/DataFrameAdapter.h"
#include "MothurDependencies/RAbundVector.h"

std::string ListVector::getOTUName(const long long bin) {
    if (static_cast<long long>(binLabels.size()) > bin) {
    } else { getLabels(); }
    return binLabels[bin];
}

ListVector::ListVector(const ListVector& other)  : DataVector(other) {
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
    Utils util;
    data.push_back(seqNames);
    const int nNames = util.getNumNames(seqNames);

    numBins++;

    if (nNames > maxRank) { maxRank = nNames; }

    numSeqs += nNames;
}

void ListVector::set(const int binNumber, const std::string &seqNames) {
    Utils util;
    const int nNames_old = util.getNumNames(data[binNumber]);
    data[binNumber] = seqNames;
    const int nNames_new = util.getNumNames(seqNames);

    if(nNames_old == 0)			{	numBins++;				}
    if(nNames_new == 0)			{	numBins--;				}
    if(nNames_new > maxRank)	{	maxRank = nNames_new;	}

    numSeqs += (nNames_new - nNames_old);
}

std::vector<std::string> ListVector::getLabels() {
    Utils util;
    util.getOTUNames(binLabels, numBins, otuTag);
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
    // output_cluster += label + "\t" + std::to_string(numBins);
    //TestHelper::Print(output_cluster);
    Utils util;
    std::vector<listCt> hold;
    for (const auto & i : data) {
        if (!i.empty()) {
            std::vector<std::string> binNames;
            std::string bin = i;
            util.splitAtComma(bin, binNames);
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
    std::sort(hold.begin(), hold.end(), abundNamesSort2); // Mothur sorts
    // its bins
    // This means that we can create an rabundvector just by sorting it by size.
    // And it should be equal. Rabund = binSize. And since it is sorted, they should be equal.

    for (auto & i : hold) {
        if (!i.bin.empty()) {
            //TestHelper::Print('\t' + hold[i].bin);
            output_cluster += "\t" + i.bin;
        }
    }
    output_cluster += "\n";
    return output_cluster;
}

int ListVector::size() {
    return static_cast<int>(data.size());
}

void ListVector::clear() {
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
            util.splitAtComma(bin, binNames);
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
    Utils util;
    for(const auto & i : data){
        const int binSize = util.getNumNames(i);
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
            Utils util;
            std::vector<listCt> hold;
            for (size_t i = 0; i < data.size(); i++) {
                if (!data[i].empty()) {
                    std::vector<std::string> binNames;
                    std::string bin = data[i];
                    util.splitAtComma(bin, binNames);
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
Rcpp::DataFrame ListVector::CreateDataFrameFromList() const {
    std::unordered_map<std::string, std::vector<std::string>> map;
    const std::vector<std::string> headers{"otu", "sequences"};
    int count = 1;
    for(const auto& bin : data) {
        if(bin.empty())
            continue;
        map[headers[0]].emplace_back("otu" + std::to_string(count++));
        map[headers[1]].emplace_back(bin);
    }
    return DataFrameAdapter::UnorderedMapToDataFrame(map);
}
