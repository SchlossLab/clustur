//
// Created by Gregory Johnson on 6/9/26.
//


#include "MothurDependencies/Utils.h"

#include "Clusters/AverageLinkage.h"
#include "Clusters/CompleteLinkage.h"
#include "Clusters/SingleLinkage.h"
#include "Clusters/WeightedLinkage.h"
#include "Clusters/Metrics/accuracy.h"
#include "Clusters/Metrics/f1score.h"
#include "Clusters/Metrics/fdr.h"
#include "Clusters/Metrics/fn.h"
#include "Clusters/Metrics/fp.h"
#include "Clusters/Metrics/fpfn.h"
#include "Clusters/Metrics/mcc.h"
#include "Clusters/Metrics/npv.h"
#include "Clusters/Metrics/ppv.h"
#include "Clusters/Metrics/sensitivity.h"
#include "Clusters/Metrics/specificity.h"
#include "Clusters/Metrics/tn.h"
#include "Clusters/Metrics/tp.h"
#include "Clusters/Metrics/tptn.h"

void Utils::mothurRandomShuffle(std::vector<int>& randomize){
    Rcpp::IntegerVector randomValues = Rcpp::wrap(randomize);
    const int size = static_cast<int>(randomize.size());
    randomValues = Rcpp::sample(randomValues, size);
    randomize = Rcpp::as<std::vector<int>>(randomValues);
}

void Utils::mothurRandomShuffle(std::vector<long long>& randomize){
    Rcpp::IntegerVector randomValues = Rcpp::wrap(randomize);
    const int size = static_cast<int>(randomize.size());
    randomValues = Rcpp::sample(randomValues, size);
    randomize = Rcpp::as<std::vector<long long>>(randomValues);
}

int Utils::getRandomIndex(const int highest){
        if (highest == 0) { return 0; }
    return static_cast<int>(R::runif(0, highest));
}
int Utils::getNumNames(std::string names){

    if(names.empty()){ return 0; }

    int count = 1;
    for_each(names.begin(), names.end(),[&count](char n){
        if(n == ','){ count++; }
    });
    return count;
}
int Utils::getOTUNames(std::vector<std::string>& currentLabels, const int numBins, const std::string& tagHeader){

        const auto currentLabelsSize = static_cast<int>(currentLabels.size());
        if (currentLabelsSize == numBins) {  return 0; }

        if (currentLabelsSize < numBins) {
            int maxLabelNumber = 0;
            const std::string snumBins = std::to_string(numBins);

            for (int i = 0; i < numBins; i++) {
                std::string binLabel = tagHeader;
                if (i < currentLabelsSize) { //label exists
                    if (getLabelTag(currentLabels[i]) == tagHeader) { //adjust 0's??
                        std::string sbinNumber = getSimpleLabel(currentLabels[i]);
                        int tempBinNumber; mothurConvert(sbinNumber, tempBinNumber);
                        if (tempBinNumber > maxLabelNumber) { maxLabelNumber = tempBinNumber; }
                        if (sbinNumber.length() < snumBins.length()) {
                            const int diff = static_cast<int>(snumBins.length() - sbinNumber.length());
                            for (int h = 0; h < diff; h++) { binLabel += "0"; }
                        }
                        binLabel += sbinNumber;
                        currentLabels[i] = binLabel;
                    }
                }else{ //create new label
                    std::string sbinNumber = std::to_string(maxLabelNumber+1); maxLabelNumber++;
                    if (sbinNumber.length() < snumBins.length()) {
                        const int diff = static_cast<int>(snumBins.length() - sbinNumber.length());
                        for (int h = 0; h < diff; h++) { binLabel += "0"; }
                    }
                    binLabel += sbinNumber;
                    currentLabels.push_back(binLabel);
                }
            }
        }
        return static_cast<int>(currentLabels.size());

}
bool Utils::mothurConvert(const std::string& item, int& num){
    if(!isNumeric1(item)) {
        return false;
    }
    num = std::stoi(item);
    return true;

}

bool Utils::mothurConvert(const std::string &item, float& num){
    if(!isNumeric1(item)) {
        return false;
    }
    num = std::stof(item);
    return true;
}
/***********************************************************************/
bool Utils::mothurConvert(const std::string &item, double& num){

    if(!isNumeric1(item)) {
        return false;
    }
    num = std::stod(item);
    return true;
}

std::string Utils::getSimpleLabel(const std::string &label){
    //remove OTU or phylo tag
        std::string newLabel1;
        for (size_t i = 0; i < label.length(); i++) {
            if(label[i]>47 && label[i]<58) { //is a digit
                newLabel1 += label[i];
            }
        }

        int num1;

        mothurConvert(newLabel1, num1);

        const std::string simple = std::to_string(num1);

        return simple;

}

std::string Utils::getLabelTag(const std::string &label){

    std::string tag;

    for (const auto n : label) {
        if(n >47 && n <58) { //is a digit
        }else {  tag += n;  }
    }

    return tag;
}

bool Utils::isNumeric1(const std::string& stringToCheck){

    bool numeric = false;

    if (stringToCheck == "") { numeric = false;  }
    else if(stringToCheck.find_first_not_of("0123456789.-") == std::string::npos) { numeric = true; }

    return numeric;

}
void Utils::splitAtComma(const std::string& s, std::vector<std::string>& container) {

        //parse string by delim and store in vector
        split(s, ',', std::back_inserter(container));

}

bool Utils::isEqual(const float num1, const float num2) {
    return std::fabs(num1-num2) <= std::fabs(static_cast<float>(num1 * 0.001));
}

float Utils::ceilDist(const float dist, const int precision){
        return static_cast<int>(ceil(dist * precision))/static_cast<float>(precision);
}

ClusterMethod * Utils::GetClusterMethod(const std::string &method, ListVector *listVector, SparseDistanceMatrix *matrix,
    RAbundVector &rAbund, const double cutoff, const double adjust) {
    if (method == "furthest")	return new CompleteLinkage(&rAbund, listVector, matrix, cutoff, method, adjust);
    if(method == "nearest") return new SingleLinkage(&rAbund, listVector, matrix, cutoff, method, adjust);
    if(method == "average")	return new AverageLinkage(&rAbund, listVector, matrix, cutoff, method, adjust);
    return new WeightedLinkage(&rAbund, listVector, matrix, cutoff, method, adjust);
}

ListVector Utils::CreateListVectorFromOtuList(const std::vector<std::string> &otuBins, const std::vector<std::string> &sequences) {
    const size_t size = otuBins.size();
    std::string lastOtu = otuBins[otuBins.size() - 1];
    lastOtu.erase(std::remove_if(lastOtu.begin(), lastOtu.end(), isalpha), lastOtu.end());
    const int otuSize = std::stoi(lastOtu) - 1;
    ListVector listVector(otuSize);
    std::string currentOtuBin = otuBins[0];
    std::string currentBinData = "";
    size_t index = 0;
    for (int i = 0; i < size; i++) {
        if (currentOtuBin == otuBins[i]) {
            currentBinData.append(sequences[i] + ",");
            continue;
        }
        currentOtuBin = otuBins[i];
        currentBinData.append(sequences[i] + ",");
        listVector.set(index++, {currentBinData.begin(), currentBinData.end() - 1});
        currentBinData = "";
    }
    return listVector;
}



void Utils::AddRowToDataFrameMap(std::unordered_map<std::string, std::vector<std::string>>& map,
                                 const std::string& data, const std::vector<std::string>& headers) {
    std::vector<std::string> splitStrings;
    splitAtComma(data, splitStrings);
    for(size_t i = 0; i < headers.size(); i++) {
        map[headers[i]].emplace_back(splitStrings[i]);
    }
}

void Utils::CheckForDistanceFileError(const std::set<std::string>& unknownNames) {
    if(!unknownNames.empty()) {
        int count = 0;
        std::string errorMessage = "These names were not found in the count table:\n";
        for(const auto& name : unknownNames) {
            errorMessage += name + "\n";
            count++;
            if(count >= 2)
                break;
        }
        Rcpp::stop(errorMessage + "Please ensure all names in the distance file are in the count table");
    }
}

ClusterMetric * Utils::GetClusterMetric(const std::string &metric) {
    if (metric == "mcc")        { return new MCC();              }
    if (metric == "sens")       { return new Sensitivity();      }
    if (metric == "spec")       { return new Specificity();      }
    if (metric == "tptn")       { return new TPTN();             }
    if (metric == "tp")         { return new TP();               }
    if (metric == "tn")         { return new TN();               }
    if (metric == "fp")         { return new FP();               }
    if (metric == "fn")         { return new FN();               }
    if (metric == "f1score")    { return new F1Score();          }
    if (metric == "accuracy")   { return new Accuracy();         }
    if (metric == "ppv")        { return new PPV();              }
    if (metric == "npv")        { return new NPV();              }
    if (metric == "fdr")        { return new FDR();              }
    if (metric == "fpfn")       { return new FPFN();             }
}

float Utils::removeConfidences(std::string& tax) {

    std::string temp = tax; float dummy; if (!hasConfidenceScore(temp, dummy)) { return 0; }

    std::string taxon;
    std::string newTax = "";
    std::string confidenceScore = "0";

    //remove last ";"
    if (tax.length() > 1) { tax = tax.substr(0, tax.length()-1); }
    std::vector<std::string> taxons; splitAtChar(tax, taxons, ';');

    for (const auto & taxonValue : taxons) {

        //if (m->getControl_pressed()) { return 0; }

        taxon = taxonValue;

        int pos = taxon.find_last_of('(');
        if (pos != -1) {
            //is it a number?
            int pos2 = taxon.find_last_of(')');
            if (pos2 != -1) {
                std::string temp = taxon.substr(pos+1, (pos2-(pos+1)));
                if (isPositiveNumeric(temp)) {
                    taxon = taxon.substr(0, pos); //rip off confidence
                    confidenceScore = temp;
                }
            }
        }
        taxon += ";";

        newTax += taxon;
    }

    tax = newTax;

    float confidence = 0; mothurConvert(confidenceScore, confidence);

    return confidence;
}

bool Utils::hasConfidenceScore(std::string& taxon, float& confidence) {
    const size_t openParen = taxon.find_last_of('(');
    const size_t closeParen = taxon.find_last_of(')');

    if ((openParen != std::string::npos) && (closeParen != std::string::npos)) {
        const std::string confidenceScore = taxon.substr(openParen+1, (closeParen-(openParen+1)));
        if (isPositiveNumeric(confidenceScore)) {  //its a confidence
            taxon = taxon.substr(0, openParen); //rip off confidence
            mothurConvert(confidenceScore, confidence);
            return true;
        }
        confidence = 0; //its part of the taxon
    } else {
        confidence = 0;
    }

    return false;

}

bool Utils::isPositiveNumeric(const std::string& stringToCheck){

    bool numeric = false;

    if (stringToCheck.empty()) { numeric = false;  }
    else if(stringToCheck.find_first_not_of("0123456789.") == std::string::npos) { numeric = true; }

    return numeric;

}

void Utils::splitAtChar(std::string& prefix, std::string& suffix, char c) {
    std::string individual = "";
    int estimLength = prefix.size();
    for(int i=0;i<estimLength;i++) {
        if(prefix[i] == c){
            suffix = prefix.substr(i+1);
            prefix = individual;
            break;
        }
        else{
            individual += prefix[i];
        }
    }
}

void Utils::splitAtChar(std::string& s, std::vector<std::string>& container, char symbol) {

    //special case to escape things
    if (symbol == '-') { splitAtDash(s, container); return; }

    //parse string by delim and store in vector
    split(s, symbol, std::back_inserter(container));
    return;

}


void Utils::splitAtDash(const std::string& estim, std::vector<std::string>& container) {
    std::string individual = "";
    int estimLength = estim.size();
    bool prevEscape = false;

    for(int i=0;i<estimLength;i++){
        if(estim[i] == '-'){
            if (prevEscape) {  individual += estim[i]; prevEscape = false;  } //add in dash because it was escaped.
            else {
                container.push_back(individual);
                individual = "";
            }
        }else if(estim[i] == '\\'){
            if (i < estimLength-1) {
                if (estim[i+1] == '-') { prevEscape=true; }  //are you a backslash before a dash, if yes ignore
                else { individual += estim[i]; prevEscape = false;  } //if no, add in
            }else { individual += estim[i]; }
        }else {
            individual += estim[i];
        }
    }

    container.push_back(individual);
}