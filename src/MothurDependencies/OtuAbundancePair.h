//
// Created by Gregory Johnson on 8/1/24.
//

#ifndef OTUABUNDANCEPAIR_H
#define OTUABUNDANCEPAIR_H
#include <string>
struct OTUAbundancePair {
    OTUAbundancePair(const int abundance, const std::string &otu, const std::string &group)
        : abundance(abundance),
          otu(otu),
          group(group) {
    }
    OTUAbundancePair() = default;

    int abundance;
    std::string otu;
    std::string group;
};
#endif //OTUABUNDANCEPAIR_H
