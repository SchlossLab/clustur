//
// Created by Gregory Johnson on 3/27/24.
//

#ifndef RSPARSEMATRIXTOCPPSPARSEMATRIXADAPTER_H
#define RSPARSEMATRIXTOCPPSPARSEMATRIXADAPTER_H
#include <vector>
#include <string>
#include <set>
#include <unordered_set>
#include <list>
#include "../MothurDependencies/OptiMatrix.h"
#include <fstream>
#include <unordered_map>
#include <sstream>


class OptimatrixAdapter {
    ///Closeness is the data represntation, it is a vector of sets that tells us which values are closes to it ->
    /// So, close[0], if Namemap[0] = A, then close[0] is the values that are within the cutoff to A.
    ///Namemap is just the names of the sequences
    ///Singletons are the values that did not make it
    ///
public:
    explicit OptimatrixAdapter(const double cutoff) {
        this->cutoff = cutoff;
    }
    OptiMatrix* ConvertToOptimatrix(const std::vector<int>
        &xPosition, const std::vector<int>& yPosition, const std::vector<double>& data);
    std::vector<std::set<long long>> GetCloseness() {return closeness;}
    std::vector<std::string> GetNameList() {return nameList;}
    std::vector<std::string> GetSingletons() {return singletons;}
private:
    double cutoff;
    //TODO Get rid of values inside of the adapter, it should not care about them
    std::vector<std::set<long long>> closeness;
    std::vector<std::string> nameList;
    std::vector<std::string> singletons;

};
#endif //RSPARSEMATRIXTOCPPSPARSEMATRIXADAPTER_H
