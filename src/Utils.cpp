//
// Created by Gregory Johnson on 3/29/24.
//

#include "MothurDependencies/Utils.h"

#include <random>
#include <sstream>
#include <unordered_set>

void Utils::mothurRandomShuffle(std::vector<int>& randomize){
    std::shuffle(randomize.begin(), randomize.end(), mersenne_twister_engine);
}
void Utils::mothurRandomShuffle(std::vector<std::string>& randomize){
    std::shuffle(randomize.begin(), randomize.end(), mersenne_twister_engine);
}
int Utils::getNumNames(std::string names){

    if(names.empty()){ return 0; }

    int count = 1;
    for_each(names.begin(), names.end(),[&count](char n){
        if(n == ','){ count++; }
    });
    return count;
}
int Utils::getOTUNames(std::vector<std::string>& currentLabels, int numBins, std::string tagHeader){

        if (currentLabels.size() == numBins) {  return 0; }

        int maxLabelNumber = 0;
        if (currentLabels.size() < numBins) {
            std::string snumBins = std::to_string(numBins);

            for (int i = 0; i < numBins; i++) {
                std::string binLabel = tagHeader;
                if (i < currentLabels.size()) { //label exists
                    if (getLabelTag(currentLabels[i]) == tagHeader) { //adjust 0's??
                        std::string sbinNumber = getSimpleLabel(currentLabels[i]);
                        int tempBinNumber; mothurConvert(sbinNumber, tempBinNumber);
                        if (tempBinNumber > maxLabelNumber) { maxLabelNumber = tempBinNumber; }
                        if (sbinNumber.length() < snumBins.length()) {
                            int diff = snumBins.length() - sbinNumber.length();
                            for (int h = 0; h < diff; h++) { binLabel += "0"; }
                        }
                        binLabel += sbinNumber;
                        currentLabels[i] = binLabel;
                    }
                }else{ //create new label
                    std::string sbinNumber = std::to_string(maxLabelNumber+1); maxLabelNumber++;
                    if (sbinNumber.length() < snumBins.length()) {
                        int diff = snumBins.length() - sbinNumber.length();
                        for (int h = 0; h < diff; h++) { binLabel += "0"; }
                    }
                    binLabel += sbinNumber;
                    currentLabels.push_back(binLabel);
                }
            }
        }
        return currentLabels.size();

}
bool Utils::mothurConvert(std::string item, int& num){
    bool error = false;

    if (isNumeric1(item)) { convert<int>(item, num); }
    else {
        num = 0;
        error = true;
    }

    return error;

}

bool Utils::mothurConvert(std::string item, float& num){
        bool error = false;

        if (isNumeric1(item)) {
            convert<float>(item, num);
        }else {
            try {
                num = std::atof(item.c_str());
            }catch(std::exception& e) {
                num = 0;
                error = true;
            }
        }

        return error;
}
/***********************************************************************/
bool Utils::mothurConvert(std::string item, double& num){

        bool error = false;

        if (isNumeric1(item)) {
            convert<double>(item, num);
        }else {
            try {
                num = std::atof(item.c_str());
            }catch(std::exception& e) {
                num = 0;
                error = true;
            }
        }

        return error;
}
template<typename T>
void Utils::convert(const std::string& s, T& x, bool failIfLeftoverChars){

    std::istringstream i(s);
    char c;
    if (!(i >> x) || (failIfLeftoverChars && i.get(c)))
        throw BadConversion(s);

}
std::string Utils::getSimpleLabel(std::string label){
        std::string simple = "";

        //remove OTU or phylo tag
        std::string newLabel1 = "";
        for (int i = 0; i < label.length(); i++) {
            if(label[i]>47 && label[i]<58) { //is a digit
                newLabel1 += label[i];
            }
        }

        int num1;

        mothurConvert(newLabel1, num1);

        simple = std::to_string(num1);

        return simple;

}

std::string Utils::getLabelTag(std::string label){

    std::string tag = "";

    for (auto n : label) {
        if(n >47 && n <58) { //is a digit
        }else {  tag += n;  }
    }

    return tag;
}

bool Utils::isNumeric1(std::string stringToCheck){

    bool numeric = false;

    if (stringToCheck == "") { numeric = false;  }
    else if(stringToCheck.find_first_not_of("0123456789.-") == std::string::npos) { numeric = true; }

    return numeric;

}
void Utils::splitAtComma(std::string& s, std::vector<std::string>& container) {

        //parse string by delim and store in vector
        split(s, ',', std::back_inserter(container));

}

void Utils::splitAtComma(std::string& prefix, std::string& suffix){
        prefix = suffix.substr(0,suffix.find_first_of(','));
        if ((suffix.find_first_of(',')+2) <= suffix.length()) {  //checks to make sure you don't have comma at end of string
            suffix = suffix.substr(suffix.find_first_of(',')+1, suffix.length());
           std::string space = " ";
            while(suffix.at(0) == ' ')
                suffix = suffix.substr(1, suffix.length());
        }else {  suffix = "";  }


}

void Utils::splitAtComma(std::string& s, std::vector<int>& container) {

    std::vector<std::string> items; splitAtComma(s, items);
    for (std::string i : items) {
        int num; mothurConvert(i, num);
        container.push_back(num);
    }

}
