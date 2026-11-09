//
// Created by Gregory Johnson on 11/6/26.
//

#include "OneGapIgnore.h"


double OneGapIgnore::Execute(const std::string &seqA, const std::string &seqB) const {

        bool overlap = false;

        const int start = setStartIgnoreTermGap(seqA, seqB, overlap);
        const int end = setEndIgnoreTermGap(seqA, seqB, overlap);
        double dist = 0;
        //non-overlapping sequences
        if (!overlap) { return 1.0000; }

        int maxMinLength = end - start;
        int difference = 0;
        bool openGapA = false;
        bool openGapB = false;

        for(int i=start;i<=end;i++){
            if(seqA[i] == '-' && seqB[i] == '-'){    maxMinLength--;    } //comparing gaps, ignore

            else if(seqB[i] != '-' && seqA[i] == '-'){ //seqB is a base, seqA is a gap
                if(!openGapA){
                    difference++;
                    openGapA = true;
                    openGapB = false;
                }else { maxMinLength--; }
            }
            else if(seqA[i] != '-' && seqB[i] == '-'){ //seqA is a base, seqB is a gap
                if(!openGapB){
                    difference++;
                    openGapA = false;
                    openGapB = true;
                }else { maxMinLength--; }
            }
            else if(seqA[i] != '-' && seqB[i] != '-'){ //both bases
                openGapA = false;
                openGapB = false;

                //no match
                if(seqA[i] != seqB[i]){ difference++; }
            }

            // dist = (double)difference / maxMinLength;
        }

    if(maxMinLength == 0)    {    dist = 1.0000;                                }
    else                    {    dist = static_cast<double>(difference) / maxMinLength;    }

    return dist;
}