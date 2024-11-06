//
// Created by Gregory Johnson on 5/1/24.
//

#include "Tests/ListVectorTestFixture.h"

bool ListVectorTestFixture::TestListVectorReturnsCorrectNumberOfBins(const std::vector<std::string>& mockListOfSequences,
    const int expectedResult) {
    Setup();
    for (const auto& sequences: mockListOfSequences) {
        listVector->push_back(sequences);
    }
    const bool result = listVector->getNumBins() == expectedResult;
    TearDown();
    return result;
}

bool ListVectorTestFixture::TestListVectorReturnsCorrectNumSeqs(const std::vector<std::string> &mockListOfSequences,
    const int expectedResult) {
    Setup();
    for (const auto& sequences: mockListOfSequences) {
        listVector->push_back(sequences);
    }
    const bool result = listVector->getNumSeqs() == expectedResult;
    TearDown();
    return result;
}

bool ListVectorTestFixture::TestListVectorReturnsCorrectMaxRank(const std::vector<std::string> &mockListOfSequences,
    const int expectedResult) {
    Setup();
    for (const auto& sequences: mockListOfSequences) {
        listVector->push_back(sequences);
    }
    const bool result = listVector->getMaxRank() == expectedResult;
    TearDown();
    return result;
}

bool ListVectorTestFixture::TestListVectorReturnsCorrectGetValue(const std::vector<std::string> &mockListOfSequences,
    const int index, const std::string &expectedResult) {
    Setup();
    for (const auto& sequences: mockListOfSequences) {
        listVector->push_back(sequences);
    }
    const bool result = listVector->get(index) == expectedResult;
    TearDown();
    return result;
}

bool ListVectorTestFixture::TestListVectorReturnsCorrectGetLabelsValue(
    const std::vector<std::string> &mockListOfSequences, const int expectedResult) {
    Setup();
    for (const auto& sequences: mockListOfSequences) {
        listVector->push_back(sequences);
    }
    const bool result = static_cast<int>(listVector->getLabels().size()) == expectedResult;
    TearDown();
    return result;
}

bool ListVectorTestFixture::TestListVectorReturnsCorrectGetOtuNames(const std::vector<std::string> &mockListOfSequences,
    const int binToTest, const std::string &expectedResult) {
    Setup();
    for (const auto& sequences: mockListOfSequences) {
        listVector->push_back(sequences);
    }
    const bool result = listVector->getOTUName(binToTest) == expectedResult;
    TearDown();
    return result;
}

bool ListVectorTestFixture::TestListVectorSetsLabelsCorrectly(const std::vector<std::string> &mockListOfLabels,
    const int expectedResult) {
    Setup();
    listVector->setLabels(mockListOfLabels);
    const bool result = static_cast<int>(listVector->getLabels().size()) == expectedResult;
    TearDown();
    return result;
}

bool ListVectorTestFixture::TestListVectorSetsPrintedLabelsCorrectly(const bool mockValue, const bool expectedResult) {
    Setup();
    const bool result = expectedResult == listVector->setPrintedLabels(mockValue);
    TearDown();
    return result;
}

bool ListVectorTestFixture::TestListVectorPushBackSetsDataCorrectly(const std::vector<std::string>& mockListOfSequences,
        const int index, const std::string& expectedResult) {
    Setup();
    for (const auto& sequences: mockListOfSequences) {
        listVector->push_back(sequences);
    }
    const bool result =listVector->get(index) == expectedResult;
    TearDown();
    return result;

}

bool ListVectorTestFixture::TestListVectorPrintDisplaysDataCorrectly(
    const std::vector<std::string> &mockListOfSequences, const std::string& expectedResult) {
    Setup();
    std::ofstream stream;
    listVector->setPrintedLabels(true);
    for (const auto& sequences: mockListOfSequences) {
        listVector->push_back(sequences);
    }
    // Get rid of ofstream in listVector
    const bool result = listVector->print(stream) == expectedResult;
    TearDown();
    return result;
}


bool ListVectorTestFixture::TestListVectorReturnsCorrectGetOtuNamesSize(const std::vector<std::string> &mockListOfSequences,
                                                                        const int binToTest, const int expectedResult) {
    Setup();
    for (const auto& sequences: mockListOfSequences) {
        listVector->push_back(sequences);
    }
    const bool result = static_cast<int>(listVector->getOTUName(binToTest).size()) == expectedResult;
    TearDown();
    return result;
}
bool ListVectorTestFixture::TestCreateDataFrameFromList(const std::string& label, const bool expectResult) {
    Setup();
    listVector->push_back("2");
    listVector->push_back("1");
    auto df = listVector->CreateDataFrameFromList();
    const std::vector<std::string> names = df.names();
    TearDown();
    return expectResult == !names.empty();
}


void ListVectorTestFixture::Setup() {
    listVector = new ListVector();
}

void ListVectorTestFixture::TearDown() {
    delete(listVector);
}
