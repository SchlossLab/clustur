//
// Created by Gregory Johnson on 6/9/26.
//

#ifndef REFACTOR_PDISTCELLMIN_H
#define REFACTOR_PDISTCELLMIN_H
struct PDistCellMin{
    unsigned long col;
    unsigned long row;
    PDistCellMin() = default;
    PDistCellMin(const unsigned long r, const unsigned long c) :  col(c), row(r) {}
};
#endif //REFACTOR_PDISTCELLMIN_H