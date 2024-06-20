//
// Created by Gregory Johnson on 3/29/24.
//

#ifndef LISTVECTOR_H
#define LISTVECTOR_H
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <map>
#include <iostream>

#include "DataVector.h"
#include "Utils.h"




class ListVector : public DataVector{
public:

    ListVector() {};
    explicit ListVector(const int n):	DataVector(), data(n, "") ,
    maxRank(0), numBins(0), numSeqs(0), otuTag("Otu"), printListHeaders(true){}
    ~ListVector(){};

    struct listCt{
        std::string bin;
        int binSize;
        std::string label;
        listCt() :  bin(""), binSize(0), label("") {};
        listCt(std::string b, int a, std::string l) :  bin(b), binSize(a), label(l) {}
    };


    static bool abundNamesSort2(listCt left, listCt right){
        if (left.bin == "") { return false; }
        if (right.bin == "") { return true; }
        if (left.binSize > right.binSize) {
            return true;
        }
        return false;
    }
    // Max rank is equal to the amount of items separated by commas
    // Number of sequences is the amount of sequences added to the list, (including those separated by delimiters)
    // Num bins is the amount of bins, or the amount of times you pushed back, the size of the list basically.
    int getNumBins()							{	return numBins;		}
    int getNumSeqs()							{	return numSeqs;		}
    int getMaxRank()							{	return maxRank;		}

    std::string get(int);
    // The amount of labels is equal to the size of the listVectorBins
    std::vector<std::string> getLabels();
    std::string getOTUName(int bin);
    void setLabels(std::vector<std::string>);
    bool setPrintedLabels(const bool pl) { printListHeaders = pl; return printListHeaders;}

    void push_back(const std::string&);
    void set(int binNumber, const std::string& seqNames);
    std::string print(std::ostream&);
    std::string print(std::ostream&, std::map<std::string, int>&);
    int size() override;
    void clear() override;
    void resize(int) override;
private:
    std::vector<std::string> data;  //data[i] is a list of names of sequences in the ith OTU.
    int maxRank = 0;
    int numBins = 0;
    int numSeqs = 0;
    std::vector<std::string> binLabels;
    std::string otuTag;
    bool printListHeaders;
    Utils util;
    void printHeaders(std::string&, std::map<std::string, int>&, bool);

public:

};



#endif //LISTVECTOR_H
