#include <iostream>
#include <fstream>
//TODO Next week: We have to separate the code out into a package! We have all the source code, now need to consider how
//it will look in a package!
#include "Adapters/OptimatrixAdapter.h"
#include "MothurDependencies/ListVector.h"
#include "TestHelpers/TestHelper.h"
#include "tests/OptiClusterTestFixture.h"
#include "Tests/OptimatrixAdapterTestFixture.h"
#if DEBUG_RCPP
#include <vector>
#include "Adapters/OptimatrixAdapter.h"
#include "MothurDependencies/OptiMatrix.h"
#include "MothurDependencies/ClusterCommand.h"

#include <Rcpp.h>
//[[Rcpp::export]]
std::string MatrixToOpiMatrixCluster(const std::vector<int> &xPosition,
                        const std::vector<int> &yPosition, const std::vector<double> &data, const double cutoff,
                        const int iterations = 2)
{
    OptimatrixAdapter adapter(cutoff);
    const auto optiMatrix = adapter.ConvertToOptimatrix(xPosition,yPosition,data);
    auto* command = new ClusterCommand();
    command->SetMaxIterations(iterations);
    return command->runOptiCluster(optiMatrix);
}
#endif
int main() {
    // UtilsTestFixture testFixture;
    // const std::string number = "123";
    // std::string exampleData = "1,2,3,4";
    // const std::vector<int> exampleIntConainer = {1,2,3,4};
    // const std::vector<std::string> exampleContainer = {"1", "2", "3", "4"};
    // bool res = testFixture.TestConvertWorksCorrectly(number, 123);
    // res = testFixture.TestIsNumericIsCorrect(number, true);
    // res = testFixture.TestMothurRandomShufflesRandomizesData(exampleIntConainer, true);
    // res = testFixture.TestMothurRandomShufflesRandomizesData({}, false);
    // res = testFixture.TestGetLabelTagsReturnsCorrectData("OTU300", "OTU");
    // res = testFixture.TestGetSimpleLabelReturnsCorrectData("OTU200", "200");
    // res = testFixture.TestMotherRandomShuffleOverloadRandomizesData(exampleContainer, true);
    // res = testFixture.TestSplitAtCommaReturnsCorrectData(exampleData, {}, exampleContainer);
    // res = testFixture.TestGetNumNamesReturnsCorrectNames("1,2,4", 3);
    // res = testFixture.TestMothurConvertStringToDoubleIsCorrect("2.00", 2, true);
    // res = testFixture.TestMothurConvertStringToFloatIsCorrect("2.00", 2, true);
    // res = testFixture.TestMothurConvertStringToIntIsCorrect("2", 2, true);


    return 0;
}