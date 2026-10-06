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
        const char sequenceOneChar = sequenceOne[i];
        const char sequenceTwoChar = sequenceTwo[i];

        if((sequenceOneChar == '-' && sequenceTwoChar == '-') || (sequenceOneChar == '.' && sequenceTwoChar == '-') ||
            (sequenceOneChar == '-' && sequenceTwoChar == '.')) {
            maxMinLength--;
            continue;
        }
        //trailing gaps, quit we already calculated all the diffs
        if(sequenceOneChar == '.' && sequenceTwoChar == '.'){ break; }
            
        if(sequenceTwoChar != '-' && (sequenceOneChar == '-' || sequenceOneChar == '.')){ //sequenceTwo is a base, sequenceOne is a gap
            if(!openGapA){
                difference++;
                openGapA = true;
                openGapB = false;
            }else { maxMinLength--; }
            continue;
        }
        if(sequenceOneChar != '-' && (sequenceTwoChar == '-' || sequenceTwoChar == '.')){ //sequenceOne is a base, sequenceTwo is a gap
            if(!openGapB){
                difference++;
                openGapA = false;
                openGapB = true;
            }else { maxMinLength--; }
            continue;
        }
        if(sequenceOneChar != '-' && sequenceTwoChar != '-'){ //both bases
            openGapA = false;
            openGapB = false;
                
            //no match
            if(sequenceOneChar != sequenceTwoChar){ difference++; }
        }
            

            
        // if (dist > cutoff) { return 1.0000; }
    }

    if(maxMinLength == 0) {
        return  1.0000;
    }
   return static_cast<double>(difference) / maxMinLength;
}
