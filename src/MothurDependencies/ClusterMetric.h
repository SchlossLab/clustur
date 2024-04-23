//
// Created by Gregory Johnson on 3/29/24.
//

#ifndef CLUSTERMETRIC_H
#define CLUSTERMETRIC_H
#include <string>
#include "Utils.h"



class ClusterMetric {
public:
    ClusterMetric(){}
    ClusterMetric(std::string n){  name = n; }
    virtual ~ClusterMetric(){};

    virtual double getValue(double, double, double, double) = 0; //tp, tn, fp, fn

    virtual std::string getName()		{	return name;        }
    virtual std::string getCitation() = 0;
    void citation() {}

protected:
    Utils util;
    std::string name;
};



#endif //CLUSTERMETRIC_H
