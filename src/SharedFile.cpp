//
// Created by Gregory Johnson on 7/31/24.
//

#include "MothurDependencies/SharedFile.h"
#include "MothurDependencies/ClusterExport.h"


Rcpp::DataFrame SharedFile::PrintData() const {
    const size_t size = tidySharedList.size();
    std::vector<std::string> groups(size);
    std::vector<std::string> otus(size);
    std::vector<double> abundanceList(size);
    int count = 0;
    for(const auto& abundances : tidySharedList) {
        groups[count] = abundances.group;
        otus[count] = abundances.OTU;
        abundanceList[count++] = abundances.groupAbundance;
    }
    return Rcpp::DataFrame::create(Rcpp::Named("samples") = groups,
                                                Rcpp::Named("otu") = otus,
                                                Rcpp::Named("abundance") = abundanceList);
}
