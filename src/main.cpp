#include <iostream>
//TODO Next week: We have to separate the code out into a package! We have all the source code, now need to consider how
//it will look in a package!
#include <vector>
#include "Adapters/OptimatrixAdapter.h"
//#include "MothurDependencies/OptiMatrix.h"
// #include "MothurDependencies/Clustercommand.h"

#include <Rcpp.h>
//[[Rcpp::export]]
std::string MatrixToOpiMatrixCluster(const std::vector<int> &xPosition,
                        const std::vector<int> &yPosition, const std::vector<double> &data, const double cutoff,
                        const int rowSize, const int colSize, const int iterations = 2)
{
    OptimatrixAdapter adapter(cutoff);
    const auto optiMatrix = adapter.ConvertToOptimatrix(xPosition,yPosition,data,rowSize,colSize);
    //auto* command = new ClusterCommand();
   // command->SetMaxIterations(iterations);
   //  command->runOptiCluster(optiMatrix);
   // return command->clusterMatrixOutput;
   return "";
}
