//
// Created by Gregory Johnson on 11/6/26.
//

#include "EachGapPairwiseDistance.h"

#include <Rcpp/exceptions.h>


double EachGapPairwiseDistance::Execute(const std::string& seqA, const std::string& seqB) const {

    int alignLength = seqA.length();
    bool overlap = false;

    int start = setStartIgnoreTermGap(seqA, seqB, overlap);
    int end = setEndIgnoreTermGap(seqA, seqB, overlap);
    double dist;

    //non-overlapping sequences
    if (!overlap) { return 1.0000; }

    int maxMinLength = end - start + 1;
    int diff = 0;
    for(int i=start;i<alignLength;i++){
        if(seqA[i] == '.' || seqB[i] == '.'){ //reached terminal gaps, so quit
            break;
        }
        else if((seqA[i] == '-' && seqB[i] == '-') || (seqA[i] == '-' && seqB[i] == '.') || (seqA[i] == '.' && seqB[i] == '-')){ maxMinLength--; } //comparing gaps, ignore
        else{
            if(seqA[i] != seqB[i]){
                diff++;
            }
        }

        // dist = (double)diff / maxMinLength;
    }

    if(maxMinLength == 0)    {    dist = 1.0000;                                }
    else            {    dist = (static_cast<double>(diff)  / static_cast<double>(maxMinLength));        }

    return dist;
}
