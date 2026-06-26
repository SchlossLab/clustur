//
// Created by Gregory Johnson on 6/9/26.
//

#ifndef REFACTOR_PDISTCELL_H
#define REFACTOR_PDISTCELL_H
struct PDistCell{
    long long index;
    float dist;
    PDistCell() :  index(0), dist(0) {};
    PDistCell(const long long c, const float d) :  index(c), dist(d) {}
    static bool CompareIndexes(const PDistCell& left, const PDistCell& right){return left.index > right.index; }
};
#endif //REFACTOR_PDISTCELL_H