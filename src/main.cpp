#include <iostream>
#include <fstream>
//TODO Next week: We have to separate the code out into a package! We have all the source code, now need to consider how
//it will look in a package!
#include "Adapters/OptimatrixAdapter.h"
#include "MothurDependencies/Cluster.h"
#include "MothurDependencies/SingleLinkage.h"
#include "TestHelpers/TestHelper.h"
#include "Tests/UtilsTestFixture.h"
#include <unordered_set>
#include "fstream"
#include "MothurDependencies/ClusterCommand.h"
#include "MothurDependencies/CompleteLinkage.h"
#include "MothurDependencies/ReadPhylipMatrix.h"
#if DEBUG_RCPP
#include <vector>
#include "Adapters/OptimatrixAdapter.h"
#include "MothurDependencies/OptiMatrix.h"
#include "MothurDependencies/ClusterCommand.h"

#include <Rcpp.h>
//[[Rcpp::export]]
std::vector<std::string> MatrixToOpiMatrixCluster(const std::vector<int> &xPosition,
                                                  const std::vector<int> &yPosition, const std::vector<double> &data,
                                                  const double cutoff,
                                                  const int iterations = 2, const bool shuffle = true) {
    OptimatrixAdapter adapter(cutoff);
    const auto optiMatrix = adapter.ConvertToOptimatrix(xPosition, yPosition, data);
    auto *command = new ClusterCommand();
    command->SetOpticlusterRandomShuffle(shuffle);
    command->SetMaxIterations(iterations);
    return command->runOptiCluster(optiMatrix);
}
#endif


std::string DistanceMatrixToDistFile(const SparseDistanceMatrix &matrix, const std::vector<std::string> &names,
                                     const std::string &outputPath) {
    int count = 0;
    const size_t size = names.size();
    std::string distanceString = "\t" + std::to_string(size) + "\n";
    for (const auto &cells: matrix.seqVec) {
        if (cells.empty()) continue;
        distanceString += names[count++];
        int indexCounter = 0;
        for (int i = 0; i < size; i++) {
            if (!cells.empty() && cells[indexCounter].index == i) {
                distanceString += "\t" + std::to_string(cells[indexCounter++].dist);
                continue;
            }
            distanceString += "\t0";
        }
        distanceString += "\n";
    }
#if !DEBUG_RCPP
        std::ofstream writeOut(outputPath);
        if (!writeOut.is_open()) {
            std::cout << "Failed to open: \n" << std::endl;
            return distanceString;
        }
        writeOut << distanceString;
        writeOut.close();
#endif
    return distanceString;
}


ReadPhylipMatrix* DistanceMatrixToRowData(const std::vector<int> &xPosition,
                                              const std::vector<int> &yPosition,
                                              const std::vector<double> &data,
                                              std::set<int> &names,
                                              const double cutoff) {
    // The indexes are +1, i need to push them back so that 1 -> 0, 2-> 1, etc (name map maybe?)

    const size_t nSeqs = data.size();
    std::unordered_map<int, RowData> dataList;
    std::unordered_map<int, int> positionsOfIndexs;
    for (size_t i = 0; i < nSeqs; i++) {
        names.insert(xPosition[i]);
        names.insert(yPosition[i]);
    }

    auto nameIter = names.begin();
    const size_t nameSize = names.size();
    for (int i = 0; i < nameSize; i++) {
        positionsOfIndexs[*nameIter] = i;
        dataList[i].name = std::to_string(*nameIter++);
        dataList[i].rowValues = std::vector<float>(nameSize);
    }

    for (size_t i = 0; i < nSeqs;  i++) {
        int xIndex = positionsOfIndexs[xPosition[i]];
        int yIndex = positionsOfIndexs[yPosition[i]];
        dataList[xIndex].rowValues[yIndex] = data[i];
        dataList[yIndex].rowValues[xIndex] = data[i];
    }
    std::vector<RowData> sequenceData(names.size());
    int index = 0;
    for (const auto &row: dataList) {
        sequenceData[index++] = row.second;
    }
    std::sort(sequenceData.begin(), sequenceData.end(), RowData::SortComparison);
    auto* matrix = new ReadPhylipMatrix("", cutoff);
    matrix->read(sequenceData);
    return matrix;
}

//[[Rcpp::export]]
std::string MatrixDataToPhylipFormat(const std::vector<int> &xPosition,
                                     const std::vector<int> &yPosition,
                                     const std::vector<double> &data,
                                     const double cutoff,
                                     const std::string &outputPath) {
    std::set<int> names;
    const auto *distanceMatrix_print = DistanceMatrixToRowData(xPosition, yPosition, data, names, cutoff);
    std::vector<std::string> nameVector(names.size());
    int count = 0;
    for (const auto &name: names) {
        nameVector[count++] = std::to_string(name);
    }
    return DistanceMatrixToDistFile(*distanceMatrix_print->getDMatrix(), nameVector, outputPath);
}

//[[Rcpp::export]]
std::string ClassicCluster(const std::vector<int> &xPosition,
                           const std::vector<int> &yPosition,
                           const std::vector<double> &data,
                           const double cutoff,
                           const std::string& method) {
    std::set<int> names;
    const auto *matrixReader = DistanceMatrixToRowData(xPosition, yPosition, data, names, cutoff);
    ClusterCommand command;
    return command.runMothurCluster(method, matrixReader->getDMatrix(), cutoff, matrixReader->getListVector());

}


void CreateRandomData(std::vector<int> &xPositions, std::vector<int> &yPositions, std::vector<double> &data,
                      const int amount) {
    xPositions.resize(amount);
    yPositions.resize(amount);
    data.resize(amount);
    std::srand(std::time(nullptr));
    for (size_t i = 0; i < amount; i++) {
        xPositions[i] = std::rand() % amount;
        yPositions[i] = i;
        const double value = (std::rand() % amount) + 1;
        data[i] = value / 250;
    }
    for (size_t i = 0; i < amount / 2; i++) {
        const int currentValue = yPositions[i];
        const int randomIndex = std::rand() % amount;
        yPositions[i] = yPositions[randomIndex];
        yPositions[randomIndex] = currentValue;
    }
}

// int main() {
//     constexpr double cutoff = 0.2;
//      std::vector<int> xPosition = {1, 2, 2, 3, 4};
//      std::vector<int> yPosition = {2, 3, 4, 4, 5};
//      std::vector<double> data = {.13f, .14f, .16f, .11f, .19f};
//
//     CreateRandonData(xPosition, yPosition,data, 100);
//     std::cout << ClassicCluster(xPosition,yPosition, data, 0.2, "furthest");
//     // ClusterClassic(xPosition,yPosition,data,cutoff);
//     // auto *distanceMatrix_print = MatrixToCluster(xPosition, yPosition, data, cutoff);
//     // OptimatrixAdapter optiMatrixAdapter(cutoff);
//     // auto *optiMatrix = optiMatrixAdapter.ConvertToOptimatrix(xPosition, yPosition, data);
//     // MatrixDataToPhylipFormat(xPosition, yPosition, data, cutoff, "/Users/grejoh/Documents/OptiClusterPackage/Opticluster/tests/output_10000.txt");
//      //istanceMatrixToRowData(xPosition, yPosition, data);
//     // const std::string path = "/Users/grejoh/Documents/OptiClusterPackage/Opticluster/tests/output.txt";
//     // ReadPhylipMatrix reader(path, cutoff);
//     // reader.read();
//     // auto *distanceMatrix = reader.getDMatrix();
//     // auto *list = reader.getListVector();
//     // auto *cluster = new SingleLinkage(new RAbundVector(), list,
//     //                                   distanceMatrix, 0.2, "nearest", -1.0);
//     // ClusterCommand command;
//     // command.runMothurCluster(cluster, distanceMatrix, cutoff, list);
//     return 0;
// }
