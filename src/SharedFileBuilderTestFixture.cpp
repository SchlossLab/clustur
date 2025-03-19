//
// Created by Gregory Johnson on 8/29/24.
//

#include "Tests/SharedFileBuilderTestFixture.h"

bool SharedFileBuilderTestFixture::TestBuildSharedFile(const ListVector& listVector,
                                                       const CountTableAdapter &countTable, const bool expectedResult) {
    Setup();
    const SharedFile* file = builder->BuildSharedFile(listVector, countTable, "otu");
    TearDown();
    return expectedResult == (file != nullptr);

}

void SharedFileBuilderTestFixture::Setup() {
    builder = new SharedFileBuilder();
}

void SharedFileBuilderTestFixture::TearDown() {
    delete(builder);
}
