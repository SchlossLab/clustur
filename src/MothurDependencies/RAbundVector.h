//
// Created by Gregory Johnson on 6/17/24.
//

#ifndef RABUNDVECTOR_H
#define RABUNDVECTOR_H



/*  Data Structure for a rabund file.
    This class is a child to datavector.  It represents OTU information at a certain distance.
    A rabundvector can be converted into and ordervector, listvector or sabundvector.
    Each member of the internal container "data" represents an individual OTU.
    So data[0] = 6, because there are six member in that OTU.
    example: listvector		=	a,b,c,d,e,f		g,h,i		j,k		l		m
             rabundvector	=	6				3			2		1		1
             sabundvector	=	2		1		1		0		0		1
             ordervector	=	1	1	1	1	1	1	2	2	2	3	3	4	5 */

#include <vector>

#include "DataVector.h"

class RAbundFloatVector;
class OrderVector;

class RAbundVector : public DataVector {

public:
    RAbundVector();
    RAbundVector(int);
    RAbundVector(std::vector<int>, int, int, int);
    //RAbundVector(vector<int>);
    RAbundVector(std::string, std::vector<int>);
    RAbundVector(const RAbundVector& bv) : DataVector(bv), data(bv.data), maxRank(bv.maxRank), numBins(bv.numBins), numSeqs(bv.numSeqs) {};
    RAbundVector(std::ifstream&);
    RAbundVector(std::ifstream& f, std::string l); //filehandle, label
    ~RAbundVector();

    int getNumBins();
    int getNumSeqs();
    int getMaxRank();

    int remove(int);
    void set(int, int);
    int get(int);
    std::vector<int> get() { return data; }
    void push_back(int);
    void pop_back();
    void resize(int);
    int size();
    void quicksort();
    int sum();
    int sum(int);
    int numNZ();
    std::vector<int> getSortedD();
    void clear();
    std::vector<int>::reverse_iterator rbegin();
    std::vector<int>::reverse_iterator rend();

    std::string print(std::ostream&); //sorted
    void nonSortedPrint(std::ostream&); //nonsorted

    RAbundVector getRAbundVector();
    RAbundFloatVector getRAbundFloatVector();
    //SAbundVector getSAbundVector();
    OrderVector getOrderVector(std::map<std::string,int>*);

private:
    std::vector<int> data;
    int maxRank;
    int numBins;
    int numSeqs;
};




#endif //RABUNDVECTOR_H
