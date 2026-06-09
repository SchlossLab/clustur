//
// Created by Gregory Johnson on 6/9/26.
//

#ifndef REFACTOR_RABUNDVECTOR_H
#define REFACTOR_RABUNDVECTOR_H




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

#include <fstream>
#include <vector>
#include <algorithm>


class RAbundVector final {

public:
    RAbundVector();
    explicit RAbundVector(int);
    ~RAbundVector();

    int getNumBins();
    int getNumSeqs();
    int getMaxRank();

    void print();

    int remove(int);
    void set(int, int);
    [[nodiscard]] int get(int) const;
    std::vector<int> get() { return data; }
    void push_back(int);
    void resize(int);
    [[nodiscard]] int size() const;
    bool quicksort();
    void clear();
    void setLabel(const std::string& l)	{label = l; }
    [[nodiscard]] std::string getLabel() const { return label; }


private:
    std::string label;
    std::vector<int> data;
    int maxRank;
    int numBins;
    int numSeqs;
};





#endif //REFACTOR_RABUNDVECTOR_H