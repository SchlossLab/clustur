//
// Created by Gregory Johnson on 11/6/26.
//

#include "MothurDependencies/IgnoreGap.h"


double IgnoreGap::Execute(const std::string &seqA, const std::string &seqB) const {

    int diff = 0;
    int start = 0;
    int end = 0;
    bool overlap = false;
    double dist = 0;

    const int alignLength = static_cast<int>(seqA.length());

    for(int i=0;i<alignLength;i++){
        if(seqA[i] != '.' && seqB[i] != '.'){
            start = i;
            overlap = true;
            break;
        }
    }

    for(int i=alignLength-1;i>=0;i--){
        if(seqA[i] != '.' && seqB[i] != '.'){
            end = i;
            overlap = true;
            break;
        }
    }

    //non-overlapping sequences
    if (!overlap) { return 1.0000; }

    int maxMinLength = end - start + 1;

    for(int i=start; i<alignLength; i++){
        if(seqA[i] == '.' || seqB[i] == '.'){
            break;
        }
        else if((seqA[i] != '-' && seqB[i] != '-')){
            if(seqA[i] != seqB[i]){
                diff++;
            }
        }else {  maxMinLength--; }

        // dist = (double)diff / maxMinLength;

    }

    if(maxMinLength == 0)    {    dist = 1.0000;                                    }
    else                    {    dist = (static_cast<double>(diff)  / static_cast<double>(maxMinLength));    }

    return dist;
}