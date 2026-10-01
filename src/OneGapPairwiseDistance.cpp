//
// Created by Gregory Johnson on 8/11/26.
//

#include "MothurDependencies/OneGapPairwiseDistance.h"
#include <Rcpp.h>

double OneGapPairwiseDistance::Execute(const std::string &sequenceOne, const std::string &sequenceTwo) const {
    int difference = 0;
    bool openGapA = false;
    bool openGapB = false;
    const int alignLength = static_cast<int>(sequenceOne.length());
        
    const int start = setStart(sequenceOne, sequenceTwo);
    const int end = setEnd(sequenceOne, sequenceTwo);
        
    int maxMinLength = end - start + 1;
        
    for(int i=start;i<alignLength;i++){
            
        //comparing gaps, ignore
        if((sequenceOne[i] == '-' && sequenceTwo[i] == '-') || (sequenceOne[i] == '.' && sequenceTwo[i] == '-') || (sequenceOne[i] == '-' && sequenceTwo[i] == '.')){    maxMinLength--;    }
        //trailing gaps, quit we already calculated all the diffs
        else if(sequenceOne[i] == '.' && sequenceTwo[i] == '.'){ break; }
            
        else if(sequenceTwo[i] != '-' && (sequenceOne[i] == '-' || sequenceOne[i] == '.')){ //sequenceTwo is a base, sequenceOne is a gap
            if(!openGapA){
                difference++;
                openGapA = true;
                openGapB = false;
            }else { maxMinLength--; }
        }
        else if(sequenceOne[i] != '-' && (sequenceTwo[i] == '-' || sequenceTwo[i] == '.')){ //sequenceOne is a base, sequenceTwo is a gap
            if(!openGapB){
                difference++;
                openGapA = false;
                openGapB = true;
            }else { maxMinLength--; }
        }
        else if(sequenceOne[i] != '-' && sequenceTwo[i] != '-'){ //both bases
            openGapA = false;
            openGapB = false;
                
            //no match
            if(sequenceOne[i] != sequenceTwo[i]){ difference++; }
        }
            

            
        // if (dist > cutoff) { return 1.0000; }
    }

    if(maxMinLength == 0) {
        return  1.0000;
    }
   return static_cast<double>(difference) / maxMinLength;
}
