//
// Created by Gregory Johnson on 3/29/24.
//
#include "MothurDependencies/ListVector.h"
#include "Adapters/DataFrameAdapter.h"
#include "MothurDependencies/RAbundVector.h"

std::string ListVector::getOTUName(const int bin) {
    if (binLabels.size() > bin) {
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

std::string ListVector::get(int index) {
    if (index < data.size()) { return data[index]; }
    return "";
}

void ListVector::setLabels(std::vector<std::string> labels) {
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
    for (int i = 0; i < data.size(); i++) {
        if (data[i] != "") {
            std::vector<std::string> binNames;
            std::string bin = data[i];
            util.splitAtComma(bin, binNames);
            int total = 0;
            for (int j = 0; j < binNames.size(); j++) {
                std::map<std::string, int>::iterator it = ct.find(binNames[j]);
                if (it == ct.end()) {
                } else { total += ct[binNames[j]]; }
            }
            listCt temp(data[i], total, "");
            hold.push_back(temp);
        }
    }
    std::sort(hold.begin(), hold.end(), abundNamesSort2); // Mothur sorts
    // its bins
    // This means that we can create an rabundvector just by sorting it by size.
    // And it should be equal. Rabund = binSize. And since it is sorted, they should be equal.

    for (int i = 0; i < hold.size(); i++) {
        if (hold[i].bin != "") {
            //TestHelper::Print('\t' + hold[i].bin);
            output_cluster += "\t" + hold[i].bin;
        }
    }
    output_cluster += "\n";
    return output_cluster;
}

int ListVector::size() {
    return data.size();
}

void ListVector::clear() {
}
void ListVector::resize(const int size) {
    data.resize(size);
}

std::string ListVector::print(std::ostream &output) {
    std::map<std::string, int> ct;
    for (int i = 0; i < data.size(); i++) {
        if (data[i] != "") {
            std::string bin = data[i];
            std::vector<std::string> binNames;
            util.splitAtComma(bin, binNames);
            if (!std::isdigit(bin[0])) {
                binNames[0] = "";
                //continue;
            }
            for (int j = 0; j < binNames.size(); j++) {
                std::string key = binNames[j];
                ct[key] = 1;
            }
        }
    }

    return print(output, ct);
}
RAbundVector ListVector::getRAbundVector() const {
    RAbundVector rav;
    Utils util;
    for(int i=0;i<data.size();i++){
        const int binSize = util.getNumNames(data[i]);
        rav.push_back(binSize);
    }
    rav.setLabel(label);

    return rav;
}

void ListVector::printHeaders(std::string &output, std::map<std::string, int> &ct, bool sortPlease) {
    if (printListHeaders) {
        if (binLabels.size() == 0) { sortPlease = false; } //we are creating arbitary otuNames
        std::vector<std::string> theseLabels = getLabels();
        if (sortPlease) {
            Utils util;
            std::vector<listCt> hold;
            for (int i = 0; i < data.size(); i++) {
                if (data[i] != "") {
                    std::vector<std::string> binNames;
                    std::string bin = data[i];
                    util.splitAtComma(bin, binNames);
                    int total = 0;
                    for (int j = 0; j < binNames.size(); j++) {
                        std::map<std::string, int>::iterator it = ct.find(binNames[j]);
                        if (it == ct.end()) {
                            //m->mothurOut("[ERROR]: " + binNames[j] + " is not in your count table. Please correct.\n"); m->setControl_pressed(true);
                        } else { total += ct[binNames[j]]; }
                    }
                    listCt temp(data[i], total, theseLabels[i]);
                    hold.push_back(temp);
                }
            }
            std::sort(hold.begin(), hold.end(), abundNamesSort2);

            //print original label for sorted by abundance otu
            for (int i = 0; i < hold.size(); i++) { output += ('\t' + hold[i].label); }
        } else {
            for (int i = 0; i < theseLabels.size(); i++) { output += ('\t' + theseLabels[i]); }
        }

        output += "\n";
        printListHeaders = false;
    }
}
///Create Test
Rcpp::DataFrame ListVector::CreateDataFrameFromList(const std::string& label) {
    std::unordered_map<std::string, std::vector<std::string>> map;
    const std::vector<std::string> headers{"otu", "bins", "label"};
    int count = 1;
    for(const auto& bin : data) {
        if(bin.empty())
            continue;
        map[headers[0]].emplace_back("otu" + std::to_string(count++));
        map[headers[1]].emplace_back(bin);
        map[headers[2]].emplace_back(label);
    }
    return DataFrameAdapter::UnorderedMapToDataFrame(map);
}