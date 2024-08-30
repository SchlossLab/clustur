//
// Created by Gregory Johnson on 8/29/24.
//

#ifndef SHAREDFILETESTFIXTURE_H
#define SHAREDFILETESTFIXTURE_H
#include "../MothurDependencies/OtuAbundancePair.h"
#include <vector>
#include <Rcpp.h>

#include "TestFixture.h"
#include "../MothurDependencies/SharedFile.h"

class SharedFileTestFixture final : public TestFixture {
public:
    bool TestSharedFilePrintData(const std::vector<SharedAbundance>& data, const Rcpp::DataFrame& expectedResult);

private:
    void Setup() override;
    void TearDown() override;
    SharedFile* sharedFile = nullptr;
};



#endif //SHAREDFILETESTFIXTURE_H
