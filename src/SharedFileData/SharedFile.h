//
// Created by Gregory Johnson on 7/31/24.
//

#ifndef SHAREDFILE_H
#define SHAREDFILE_H
#include <vector>
#include <Rcpp.h>
#include "../DataStructures/OtuAbundancePair.h"


class SharedFile {
public:
    explicit SharedFile(const std::vector<SharedAbundance>& otherTidySharedList)
        : tidySharedList(otherTidySharedList) {
    }
    Rcpp::DataFrame PrintData(const std::string &binName) const;
private:
    std::vector<SharedAbundance> tidySharedList;
};



#endif //SHAREDFILE_H
