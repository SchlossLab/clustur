//
// Created by Gregory Johnson on 11/6/26.
//

#include "MothurDependencies/EachGapPairwiseDistance.h"


double EachGapPairwiseDistance::Execute(const std::string& seqA, const std::string& seqB,
    const double cutoff) const {

    const int alignLength = static_cast<int>(seqA.length());

    const int start = setStart(seqA, seqB);
    const int end = setEnd(seqA, seqB);

    int maxMinLength = end - start + 1;
    int diff = 0;
    double dist = 0;

    for(int i=start;i<alignLength;i++){
        if(seqA[i] == '.' && seqB[i] == '.'){ //reached terminal gaps, so quit
            break;
        }
        else if((seqA[i] == '-' && seqB[i] == '-') || (seqA[i] == '-' && seqB[i] == '.') || (seqA[i] == '.' && seqB[i] == '-')){ maxMinLength--; } //comparing gaps, ignore
        else{
            if(seqA[i] != seqB[i]){
                diff++;
            }
        }

        dist = static_cast<double>(diff) / maxMinLength;
        if (dist > cutoff) return 1.0000;
    }

    if(maxMinLength == 0)  return  1.0000;
    return (static_cast<double>(diff)  / static_cast<double>(maxMinLength));
}
