//
// Created by Gregory Johnson on 3/29/24.
//
#include "MothurDependencies/ListVector.h"

std::string ListVector::getOTUName(int bin) {
    if (binLabels.size() > bin) {
    } else { getLabels(); }
    return binLabels[bin];
}


/***********************************************************************/

void ListVector::push_back(std::string seqNames) {
    Utils util;
    data.push_back(seqNames);
    int nNames = util.getNumNames(seqNames);

    numBins++;

    if (nNames > maxRank) { maxRank = nNames; }

    numSeqs += nNames;
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
    printListHeaders = true;
    otuTag = "Otu";
    std::string output_cluster;
    printHeaders(output_cluster, ct, true);
    //output_cluster += label + "\t" + std::to_string(numBins);
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
    std::sort(hold.begin(), hold.end(), abundNamesSort2);

    for (int i = 0; i < hold.size(); i++) {
        if (hold[i].bin != "") {
            //TestHelper::Print('\t' + hold[i].bin);
            output_cluster += "\t" + hold[i].bin;
        }
    }
    output_cluster += "\n";
    return output_cluster;
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
