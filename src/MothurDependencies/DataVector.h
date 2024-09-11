//
// Created by Gregory Johnson on 6/17/24.
//

#ifndef DATAVECTOR_H
#define DATAVECTOR_H
#include "Utils.h"
#include <map>
#include <fstream>



class DataVector {

public:
    virtual ~DataVector() = default;

    DataVector(){ }
    DataVector(const std::string& l) : label(l) {};
    virtual int size() = 0;
    virtual void clear() = 0;

    virtual void resize(int) = 0;

    virtual std::string print(std::ostream&, std::map<std::string, int>&) {return "";}
    virtual void printHeaders(std::string&, std::map<std::string, int>&, bool) {};
    virtual std::string print(std::ostream&, bool&) {return "";}
    virtual std::string print(std::ostream&) {return "";}
    virtual void print() {};

    void setLabel(const std::string& l)		{	label = l;			}
    std::string getLabel()			    {	return label;		}


protected:
    std::string label;
    Utils util;

};

#endif //DATAVECTOR_H
