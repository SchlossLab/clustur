//
// Created by Gregory Johnson on 9/11/26.
//
#include "DataStructures/ClusterParameters.h"

void ClusterParameters::SetParameters(const std::string &key, const std::string &value) {
    parameters[key] = value;
}
