//
// Created by Gregory Johnson on 7/31/24.
//

#ifndef SHAREDFILEBUILDER_H
#define SHAREDFILEBUILDER_H

#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <set>
#include <list>

#include "ClusterExport.h"
#include "ListVector.h"
#include "RSparseMatrix.h"
#include "SharedFile.h"
#include "OtuAbundancePair.h"

class SharedFileBuilder {
public:
    SharedFile *BuildSharedFile(RAbundVector *rAdbundVector, const ClusterExport *clusterExport);
    SharedFileBuilder() = default;
};





#endif //SHAREDFILEBUILDER_H
