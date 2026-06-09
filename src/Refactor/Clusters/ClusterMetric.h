//
// Created by Gregory Johnson on 3/29/24.
//

#ifndef CLUSTERMETRIC_H
#define CLUSTERMETRIC_H
#include <string>
#include "../MothurDependencies/Utils.h"



class ClusterMetric {
public:
    ClusterMetric() = default;
    explicit ClusterMetric(const std::string &n){name = n; }
    virtual ~ClusterMetric() = default;
    virtual double getValue(double, double, double, double) = 0; //tp, tn, fp, fn
    virtual std::string getName() {return name;}
    virtual std::string getCitation() = 0;

protected:
    Utils util;
    std::string name;
};



#endif //CLUSTERMETRIC_H
