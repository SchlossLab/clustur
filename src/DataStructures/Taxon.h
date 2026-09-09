//
// Created by Gregory Johnson on 9/3/26.
//

#ifndef REFACTOR_TAXON_H
#define REFACTOR_TAXON_H
#include <string>
struct Taxon {
    std::string name;
    float confidence;

    Taxon(std::string n, const float conf) : name(std::move(n)), confidence(conf) {}
    ~Taxon()= default;
};
#endif //REFACTOR_TAXON_H