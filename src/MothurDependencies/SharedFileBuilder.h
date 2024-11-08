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
#include "SharedFile.h"
#include "OtuAbundancePair.h"
#include "../Adapters/CountTableAdapter.h"


class SharedFileBuilder {
public:
    SharedFile *BuildSharedFile(const ListVector &listVector,
        const CountTableAdapter& countTable);
    SharedFileBuilder() = default;
private:
    struct SampleInformation {
        SampleInformation() = default;
        SampleInformation(std::string otu, std::string group, const double abundance)
            : otu(std::move(otu)),
              group(std::move(group)),
              abundance(abundance) {
        }
        std::string otu;
        std::string group;
        double abundance = 0;
    };
};





#endif //SHAREDFILEBUILDER_H
