#include <iostream>
#include <fstream>
//TODO Next week: We have to separate the code out into a package! We have all the source code, now need to consider how
//it will look in a package!
#include "Adapters/OptimatrixAdapter.h"
#include "MothurDependencies/Cluster.h"
#include "MothurDependencies/SingleLinkage.h"
#include "TestHelpers/TestHelper.h"
#include "Tests/UtilsTestFixture.h"
#include "fstream"
#include "MothurDependencies/ClusterCommand.h"
#include "MothurDependencies/ReadPhylipMatrix.h"
#if DEBUG_RCPP
#include <vector>
#include "Adapters/OptimatrixAdapter.h"
#include "MothurDependencies/OptiMatrix.h"
#include "MothurDependencies/ClusterCommand.h"

#include <Rcpp.h>
//[[Rcpp::export]]
std::vector<std::string> MatrixToOpiMatrixCluster(const std::vector<int> &xPosition,
                        const std::vector<int> &yPosition, const std::vector<double> &data, const double cutoff,
                        const int iterations = 2, const bool shuffle = true)
{
    OptimatrixAdapter adapter(cutoff);
    const auto optiMatrix = adapter.ConvertToOptimatrix(xPosition,yPosition,data);
    auto* command = new ClusterCommand();
    command ->SetOpticlusterRandomShuffle(shuffle);
    command->SetMaxIterations(iterations);
    return command->runOptiCluster(optiMatrix);
}
#endif

SparseDistanceMatrix *MatrixToCluster(const std::vector<int> &xPosition,
                                      const std::vector<int> &yPosition, const std::vector<double> &data,
                                      const double cutoff) {
    Rcpp::Rcout << "Made it 1" << std::endl;
    const auto count = std::max(*std::max_element(xPosition.begin(), xPosition.end()),
                                *std::max_element(yPosition.begin(), yPosition.end()));
    Rcpp::Rcout << "Made it 2" << std::endl;
    std::vector<PDistCell> cells(count + 1);
    auto *distanceMatrix = new SparseDistanceMatrix();
    distanceMatrix->resize(count);
    for (size_t i = 0; i < data.size(); i++) {
        try {
            cells[i] = PDistCell(xPosition[i], data[i]);
            distanceMatrix->addCell(xPosition[i], PDistCell(yPosition[i], data[i]));
        }
        catch(std::exception& ex){
            Rcpp::Rcout << ex.what() << std::endl;
            Rcpp::Rcout << "at i postion: " << i << " : Count: " << count << std::endl;
        }
        //The rows in the sparse Distance matrix is the I'positons and the index in the Pdist cell is the j postion
        //This is why we are having duplicates and other isuses...
    }
    Rcpp::Rcout << "Made it 3" << std::endl;
    return distanceMatrix;
}

std::string DistanceMatrixToDistFile(const SparseDistanceMatrix &matrix, const std::vector<std::string> &names,
                              const std::string &outputPath) {
    int count = 0;
    const size_t size = names.size();
    Rcpp::Rcout << "Size: " << size << std::endl;
    std::string distanceString = "\t" + std::to_string(size) + "\n";
    for (const auto &cells: matrix.seqVec) {
        distanceString += names[count++];
        int indexCounter = 0;
        for (int i = 0; i < size; i++) {
            if (!cells.empty() && cells[indexCounter].index == i) {
                Rcpp::Rcout << "Adding too: "<< std::endl;
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

//[[Rcpp::export]]
std::string MatrixDataToPhylipFormat(const std::vector<int> &xPosition,
                              const std::vector<int> &yPosition,
                              const std::vector<double> &data,
                              const double cutoff,
                              const std::string &outputPath) {
    Rcpp::Rcout << "Failing here" << std::endl;
    const auto *distanceMatrix_print = MatrixToCluster(xPosition, yPosition, data, cutoff);
    OptimatrixAdapter optiMatrixAdapter(cutoff);
    Rcpp::Rcout << "Failing here1" << std::endl;
    auto *optiMatrix = optiMatrixAdapter.ConvertToOptimatrix(xPosition, yPosition, data);
    Rcpp::Rcout << "Failing here2" << std::endl;
    return DistanceMatrixToDistFile(*distanceMatrix_print, optiMatrix->GetNameList(), outputPath);
}

//
// int main() {
//     double cutoff = 0.2;
//     const std::vector<int> xPosition = {1, 2, 2, 3, 4};
//     const std::vector<int> yPosition = {2, 3, 4, 4, 5};
//     const std::vector<double> data = {.13f, .24f, .41f, .55f, .19f};
//     auto *distanceMatrix_print = MatrixToCluster(xPosition, yPosition, data, cutoff);
//     OptimatrixAdapter optiMatrixAdapter(cutoff);
//     auto *optiMatrix = optiMatrixAdapter.ConvertToOptimatrix(xPosition, yPosition, data);
//     MatrixDataToPhylipFormat(xPosition, yPosition, data, cutoff, "/Users/grejoh/Documents/OptiClusterPackage/Opticluster/tests/output_3.txt");
//
//
//     // const std::string path = "/Users/grejoh/Documents/OptiClusterPackage/Opticluster/tests/output_2.txt";
//     // ReadPhylipMatrix reader(path, cutoff);
//     // reader.read();
//     // auto *distanceMatrix = reader.getDMatrix();
//     // auto *list = reader.getListVector();
//     // auto *cluster = new SingleLinkage(new RAbundVector(), list,
//     //                                   distanceMatrix, 0.2, "nearest", -1.0);
//     // ClusterCommand command;
//     // command.runMothurCluster(cluster, distanceMatrix, cutoff, list);
//     // return 0;
// }
