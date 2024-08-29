//
// Created by Gregory Johnson on 3/29/24.
//

#ifndef UTILS_H
#define UTILS_H
#define MOTHURMAX 1e6

#include <vector>
#include <random>
#include <algorithm>
#include <set>
#include <sstream>
#include <stdexcept>

#include "PDistCell.h"

class BadConversion : public std::runtime_error {
public:
    BadConversion(const std::string& s) : runtime_error(s){ }
};

class Utils {
public:
    void mothurRandomShuffle(std::vector<int>& randomize);
    void mothurRandomShuffle(std::vector<std::string>& randomize);
    void mothurRandomShuffle(std::vector<PDistCellMin>& randomize);
    int getRandomIndex(int highest);
    int getNumNames(std::string names);
    bool mothurConvert(std::string item, int& num);
    bool isNumeric1(std::string stringToCheck);
    int getOTUNames(std::vector<std::string>& currentLabels, int numBins, std::string tagHeader);
    bool mothurConvert(std::string item, double& num);
    bool mothurConvert(std::string item, float& num);
    std::string getLabelTag(std::string label);
    std::string getSimpleLabel(std::string label);
    void splitAtComma(std::string& s, std::vector<std::string>& container);
    bool isEqual(float num1, float num2);

    float ceilDist(float dist, int precision);

    void AddRowToDataFrameMap(std::unordered_map<std::string, std::vector<std::string>> &map, std::string &data,
                              const std::vector<std::string> &headers);

    template <typename Out>
    void split(const std::string &s, char delim, Out result) {
        std::istringstream iss(s);
        std::string item;
        while (getline(iss, item, delim)) {
            if (!item.empty()) { //ignore white space
                *result++ = item;
            }
        }
    }

    //template<typename T>
    // void convert(const std::string& s, T& x, bool failIfLeftoverChars = true);
private:
    std::mt19937_64 mersenne_twister_engine;
};



#endif //UTILS_H
