//
// Created by Gregory Johnson on 9/10/26.
//

#ifndef REFACTOR_OPTIDATACOMPONENTS_H
#define REFACTOR_OPTIDATACOMPONENTS_H
#include "SparseDistanceMatrix.h"
#include "ListVector.h"
struct OptiDataComponent {
    SparseDistanceMatrix* matrix;
    ListVector* listVector;
};
#endif //REFACTOR_OPTIDATACOMPONENTS_H