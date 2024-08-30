//
// Created by Gregory Johnson on 7/31/24.
//

#ifndef SHAREDFILEBUILDER_H
#define SHAREDFILEBUILDER_H

#include <utility>
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
#include "../Adapters/CountTableAdapter.h"


class SharedFileBuilder {
public:
    SharedFile *BuildSharedFile(const std::unordered_map<std::string, ListVector> &listVectorMap,
        const CountTableAdapter& countTable);
    SharedFileBuilder() = default;
};





#endif //SHAREDFILEBUILDER_H
