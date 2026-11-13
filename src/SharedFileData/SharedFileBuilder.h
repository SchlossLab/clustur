//
// Created by Gregory Johnson on 7/31/24.
//

#ifndef SHAREDFILEBUILDER_H
#define SHAREDFILEBUILDER_H

#include <utility>
#include <vector>
#include <fstream>
#include "../Adapters/CountTableAdapter.h"
#include "../DataStructures/ListVector.h"
#include "SharedFile.h"

class SharedFileBuilder {
public:
    static SharedFile *BuildSharedFile(const ListVector &listVector,
                                       const CountTableAdapter& countTable, const std::string &binName);
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
