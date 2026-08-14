//
// Created by Gregory Johnson on 8/11/26.
//

#include "MothurDependencies/OneGapPairwiseDistance.h"
#include <Rcpp.h>

double OneGapPairwiseDistance::Execute(const std::string &sequenceOne, const std::string &sequenceTwo) {
    if (sequenceOne.size() != sequenceTwo.size()) Rcpp::stop("Two sequences must have same size");

    bool isOpenGap = false;
    size_t seqWithoutGapsSize = 0;
    int mismatch = 0;
    for (size_t i = 0; i < sequenceOne.size(); i++) {
        if (sequenceOne[i] == '-' || sequenceTwo[i] == '-') {
            if (!isOpenGap) {
                mismatch++;
                seqWithoutGapsSize++;
            }
            isOpenGap = true;
            continue;
        }
        isOpenGap = false;
        seqWithoutGapsSize++;
        if (sequenceOne[i] != sequenceTwo[i]) mismatch++;
    }
    return static_cast<double>(mismatch)/static_cast<double>(seqWithoutGapsSize);
}
