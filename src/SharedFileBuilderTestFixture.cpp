//
// Created by Gregory Johnson on 8/29/24.
//

#include "Tests/SharedFileBuilderTestFixture.h"

bool SharedFileBuilderTestFixture::TestBuildSharedFile(const std::unordered_map<std::string, ListVector>& map,
                                                       const CountTableAdapter &countTable, const bool expectedResult) {
    Setup();
    const SharedFile* file = builder->BuildSharedFile(map, countTable);
    TearDown();
    return expectedResult == (file != nullptr);

}

void SharedFileBuilderTestFixture::Setup() {
    builder = new SharedFileBuilder();
}

void SharedFileBuilderTestFixture::TearDown() {
    delete(builder);
}
