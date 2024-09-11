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
#include <fstream>
#include <vector>
#include <algorithm>

class RAbundFloatVector;
class OrderVector;

class RAbundVector final : public DataVector {

public:
    RAbundVector();
    explicit RAbundVector(int);
    ~RAbundVector() override;

    int getNumBins();
    int getNumSeqs();
    int getMaxRank();

    void print() override;

    int remove(int);
    void set(int, int);
    int get(int) const;
    std::vector<int> get() { return data; }
    void push_back(int);
    void resize(int) override;
    int size() override;
    bool quicksort();
    void clear() override;

private:
    std::vector<int> data;
    int maxRank;
    int numBins;
    int numSeqs;
};




#endif //RABUNDVECTOR_H
