//
// Created by Gregory Johnson on 6/9/26.
//

#ifndef REFACTOR_UTILS_H
#define REFACTOR_UTILS_H
#include <vector>
#include <random>
#include <set>
#include <sstream>
#include <unordered_map>
#include <Rcpp.h>
#include "../Clusters/ClusterMetric.h"

class ClusterMethod;
class SparseDistanceMatrix;
class ListVector;
class RAbundVector;
class Utils {
public:
    Utils() = default;
    static void mothurRandomShuffle(std::vector<int>& randomize);
    static void mothurRandomShuffle(std::vector<long long>& randomize);
    static int getRandomIndex(int highest);
    static int getNumNames(std::string names);
    static bool mothurConvert(const std::string& item, int& num);
    static bool isNumeric1(const std::string& stringToCheck);
    static int getOTUNames(std::vector<std::string>& currentLabels, int numBins, const std::string& tagHeader);
    static bool mothurConvert(const std::string &item, double& num);
    static bool mothurConvert(const std::string &item, float& num);
    static std::string getLabelTag(const std::string &label);
    static std::string getSimpleLabel(const std::string &label);
    static void splitAtComma(const std::string& s, std::vector<std::string>& container);
    static bool isEqual(float num1, float num2);
    static float ceilDist(float dist, int precision);
    static ClusterMethod* GetClusterMethod(const std::string& method, ListVector* listVector,
    SparseDistanceMatrix* matrix, RAbundVector& rAbund, double cutoff, double adjust = -1);
    template <typename T>
    static void SortVector(std::vector<T>&);

    static ListVector CreateListVectorFromOtuList(const std::vector<std::string> &otuBins, const std::vector<std::string> &sequences);
    static void AddRowToDataFrameMap(std::unordered_map<std::string, std::vector<std::string>> &map, const std::string &data,
                              const std::vector<std::string> &headers);
    //Error Checking
    static void CheckForDistanceFileError(const std::set<std::string>&);
    static ClusterMetric* GetClusterMetric(const std::string& metric);

    template <typename Out>
    static void split(const std::string &s, const char delim, Out result) {
        std::istringstream iss(s);
        std::string item;
        while (getline(iss, item, delim)) {
            if (!item.empty()) { //ignore white space
                *result++ = item;
            }
        }
    }
};

template<typename T>
void Utils::SortVector(std::vector<T>& vector) {
    std::sort(vector.begin(), vector.end());
}


#endif //REFACTOR_UTILS_H
