//
// Created by Gregory Johnson on 3/29/24.
//

#include "MothurDependencies/Utils.h"

#include <cmath>
#include <random>
#include <sstream>
#include <unordered_set>


void Utils::mothurRandomShuffle(std::vector<int>& randomize){
    Rcpp::IntegerVector randomValues = Rcpp::wrap(randomize);
    const int size = static_cast<int>(randomize.size());
    randomValues = Rcpp::sample(randomValues, size);
    randomize = Rcpp::as<std::vector<int>>(randomValues);
}
void Utils::mothurRandomShuffle(std::vector<std::string>& randomize){
    Rcpp::CharacterVector randomValues = Rcpp::wrap(randomize);
    const int size = static_cast<int>(randomize.size());
    randomValues = Rcpp::sample(randomValues, size);
    randomize = Rcpp::as<std::vector<std::string>>(randomValues);
}

void Utils::mothurRandomShuffle(std::vector<PDistCellMin> &randomize) {
    std::shuffle (randomize.begin(), randomize.end(), mersenne_twister_engine);
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

void Utils::AddRowToDataFrameMap(std::unordered_map<std::string, std::vector<std::string>>& map,
                                 const std::string& data, const std::vector<std::string>& headers) {
    Utils utils;
    std::vector<std::string> splitStrings;
    utils.splitAtComma(data, splitStrings);
    for(size_t i = 0; i < headers.size(); i++) {
        map[headers[i]].emplace_back(splitStrings[i]);
    }
}

void Utils::CheckForDistanceFileError(const std::set<std::string>& unknownNames) const {
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