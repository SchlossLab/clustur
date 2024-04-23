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
#include "Utils.h"




class ListVector {
public:

    ListVector() = default;
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
    int getNumBins()							{	return numBins;		}
    int getNumSeqs()							{	return numSeqs;		}
    int getMaxRank()							{	return maxRank;		}

    std::string get(int);
    std::vector<std::string> getLabels();
    std::string getOTUName(int bin);
    void setLabels(std::vector<std::string>);

    void setPrintedLabels(bool pl) { printListHeaders = pl; }

    void push_back(std::string);

    std::string print(std::ostream&);
    std::string print(std::ostream&, std::map<std::string, int>&);
    void setLabel(std::string l)		{	label = l;			}
private:
    std::vector<std::string> data;  //data[i] is a list of names of sequences in the ith OTU.
    int maxRank;
    int numBins;
    int numSeqs;
    std::vector<std::string> binLabels;
    std::string otuTag;
    bool printListHeaders;
    Utils util;
    std::string label;
    void printHeaders(std::string&, std::map<std::string, int>&, bool);

};



#endif //LISTVECTOR_H
