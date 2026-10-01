//
// Created by Gregory Johnson on 8/11/26.
//

#ifndef REFACTOR_PAIRWISEDISTANCECALCULATOR_H
#define REFACTOR_PAIRWISEDISTANCECALCULATOR_H
#include <string>
#include <vector>

class PairwiseDistanceCalculator {
public:
    PairwiseDistanceCalculator() = default;

    virtual ~PairwiseDistanceCalculator() = default;
    virtual double Execute(const std::string& sequenceOne, const std::string& sequenceTwo) const = 0;
protected:
    [[nodiscard]] int setStart(const std::string &seqA, const std::string &seqB) const {
        int start = 0;
        const int alignLength = static_cast<int>(seqA.length());
        for(int i=0;i<alignLength;i++){
            if((seqA[i] != '.' || seqB[i] != '.')){ //one of you is not a terminal gap
                start = i;
                break;
            }
        }
        return start;
    }
    /***********************************************************************/
    [[nodiscard]] int setEnd(const std::string &seqA, const std::string &seqB) const {
        int end = 0;
        const int alignLength = static_cast<int>(seqA.length());

        for(int i=alignLength-1;i>=0;i--){
            if((seqA[i] != '.' || seqB[i] != '.')){ //one of you is not a terminal gap
                end = i;
                break;
            }
        }
        return end;
    }
};


#endif //REFACTOR_PAIRWISEDISTANCECALCULATOR_H