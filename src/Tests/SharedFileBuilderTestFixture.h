//
// Created by Gregory Johnson on 8/29/24.
//

#ifndef SHAREDFILEBUILDERTESTFIXTURE_H
#define SHAREDFILEBUILDERTESTFIXTURE_H
#include "../MothurDependencies/OtuAbundancePair.h"
#include <Rcpp.h>

#include "TestFixture.h"
#include "../MothurDependencies/SharedFileBuilder.h"

class SharedFileBuilderTestFixture final : public TestFixture {
public:
    //Test if null, we are going to add error checking
    bool TestBuildSharedFile(const std::unordered_map<std::string, ListVector>&, const CountTableAdapter& countTable,
         bool expectedResult);
private:
    void Setup() override;
    void TearDown() override;
    SharedFileBuilder* builder = nullptr;
};



#endif //SHAREDFILEBUILDERTESTFIXTURE_H
