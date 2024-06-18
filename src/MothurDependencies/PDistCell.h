//
// Created by Gregory Johnson on 6/17/24.
//

#ifndef PDISTCELL_H
#define PDISTCELL_H

struct PDistCell{
    unsigned long index;
    float dist;
    PDistCell() :  index(0), dist(0) {};
    PDistCell(const unsigned long c, const float d) :  index(c), dist(d) {}
    static bool CompareIndexes(const PDistCell& left, const PDistCell& right){return left.index > right.index; }
};
struct PDistCellMin{
    unsigned long row;
    unsigned long col;

    PDistCellMin(const unsigned long r, const unsigned long c) :  col(c), row(r) {}
};
#endif //PDISTCELL_H
