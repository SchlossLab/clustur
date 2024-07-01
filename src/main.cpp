#include <iostream>
#include <fstream>
//TODO Next week: We have to separate the code out into a package! We have all the source code, now need to consider how
//it will look in a package!
#include "Adapters/OptimatrixAdapter.h"
#include "TestHelpers/TestHelper.h"
#include "Tests/UtilsTestFixture.h"
#include <unordered_set>
#include "Adapters/MatrixAdapter.h"
#include "MothurDependencies/ClusterCommand.h"
#if DEBUG_RCPP
#include <vector>
#include "Adapters/OptimatrixAdapter.h"
#include "MothurDependencies/OptiMatrix.h"
#include "MothurDependencies/ClusterCommand.h"

#include <Rcpp.h>

//[[Rcpp::export]]
void WriterPhylipFile(const std::vector<int> &xPosition,
                      const std::vector<int> &yPosition, const std::vector<double> &data,
                      const double cutoff, const std::string& saveLocation) {
    MatrixAdapter adapter(xPosition, yPosition, data, cutoff);
    adapter.CreatePhylipFile(saveLocation);
}

//[[Rcpp::export]]
std::vector<std::string> MatrixToOpiMatrixCluster(const std::vector<int> &xPosition,
                                                  const std::vector<int> &yPosition, const std::vector<double> &data,
                                                  const double cutoff,
                                                  const int maxIterations = 100, const bool shuffle = true) {
    OptimatrixAdapter adapter(cutoff);
    const auto optiMatrix = adapter.ConvertToOptimatrix(xPosition, yPosition, data);
    auto *command = new ClusterCommand();
    command->SetOpticlusterRandomShuffle(shuffle);
    command->SetMaxIterations(maxIterations);
    return command->runOptiCluster(optiMatrix);
}
#endif

//[[Rcpp::export]]
std::string ClassicCluster(const std::vector<int> &xPosition,
                           const std::vector<int> &yPosition,
                           const std::vector<double> &data,
                           const double cutoff,
                           const std::string& method) {
    MatrixAdapter adapter(xPosition, yPosition, data, cutoff);
    ClusterCommand command;

    return command.runMothurCluster(method, adapter.CreateSparseMatrix(), cutoff, adapter.GetListVector());

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
        const double value = (std::rand() % amount);
        data[i] = value / amount;
    }
    for (size_t i = 0; i < amount / 2; i++) {
        const int currentValue = yPositions[i];
        const int randomIndex = std::rand() % amount;
        yPositions[i] = yPositions[randomIndex];
        yPositions[randomIndex] = currentValue;
    }
}

int main() {
    const auto start = std::chrono::high_resolution_clock::now();
    constexpr double cutoff = 0.5;
     std::vector<int> xPosition = {1, 2, 2, 3, 4};
     std::vector<int> yPosition = {2, 3, 4, 4, 5};
     std::vector<double> data = {.13f, .14f, .16f, .11f, .19f};

    CreateRandomData(xPosition, yPosition,data, 1000);
    std::cout << ClassicCluster(xPosition,yPosition, data, cutoff, "furthest") << std::endl;
    const auto stop = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken: " << std::chrono::duration_cast<chrono::milliseconds>(stop - start).count();
    return 0;
}
//TODO: Get smallest cell and mothurRandomShuffle need to be improved
