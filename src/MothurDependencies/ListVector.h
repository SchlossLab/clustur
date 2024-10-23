//
// Created by Gregory Johnson on 3/29/24.
//

#ifndef LISTVECTOR_H
#define LISTVECTOR_H
#include <string>
#include <utility>
#include <vector>
#include <algorithm>
#include <fstream>
#include <map>
#include <iostream>
#include <Rcpp.h>
#include "DataVector.h"
#include "RAbundVector.h"
#include "Utils.h"




class ListVector final : public DataVector{
public:

    ListVector() = default;
    ListVector(const ListVector& other);
    explicit ListVector(const int n):	DataVector(), data(n, "") ,
    maxRank(0), numBins(0), numSeqs(0), otuTag("Otu"), printListHeaders(true){}
    ~ListVector() override = default;

    struct listCt{
        std::string bin;
        int binSize;
        std::string label;
        listCt() : binSize(0) {};
        listCt(std::string b, const int a, std::string l) :  bin(std::move(b)), binSize(a), label(std::move(l)) {}
    };


    static bool abundNamesSort2(const listCt& left, const listCt& right){
        if (left.bin.empty()) { return false; }
        if (right.bin.empty()) { return true; }
        if (left.binSize > right.binSize) {
            return true;
        }
        return false;
    }
    // Max rank is equal to the amount of items separated by commas
    // Number of sequences is the amount of sequences added to the list, (including those separated by delimiters)
    // Num bins is the amount of bins, or the amount of times you pushed back, the size of the list basically.
    int getNumBins() const							{	return numBins;		}
    int getNumSeqs() const							{	return numSeqs;		}
    int getMaxRank() const							{	return maxRank;		}

    std::string get(long long) const;
    // The amount of labels is equal to the size of the listVectorBins
    std::vector<std::string> getLabels();
    std::string getOTUName(long long bin);
    void setLabels(const std::vector<std::string> &);
    bool setPrintedLabels(const bool pl) { printListHeaders = pl; return printListHeaders;}

    void push_back(const std::string&);
    void set(int binNumber, const std::string& seqNames);
    std::string print(std::ostream&) override;

    RAbundVector getRAbundVector() const;

    std::string print(std::ostream&, std::map<std::string, int>&) override;
    int size() override;
    void clear() override;
    void resize(int) override;
    Rcpp::DataFrame CreateDataFrameFromList() const;

private:
    std::vector<std::string> data;  //data[i] is a list of names of sequences in the ith OTU.
    int maxRank = 0;
    int numBins = 0;
    int numSeqs = 0;
    std::vector<std::string> binLabels;
    std::string otuTag;
    bool printListHeaders{};
    Utils util;
    void printHeaders(std::string&, std::map<std::string, int>&, bool) override;



public:

};



#endif //LISTVECTOR_H
