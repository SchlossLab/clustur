#include <testthat.h>
#include "Tests/OptimatrixAdapterTestFixture.h"
// Normally this would be a function from your package's
// compiled library -- you might instead just include a header
// file providing the definition, and let R CMD INSTALL
// handle building and linking.

// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.
//  OptiMatrix* ConvertToOptimatrix(const std::vector<int>
//        &xPosition, const std::vector<int>& yPosition, const std::vector<double>& data, int rowSize, int colSize);
void Check()
{
   
}
context("Optimatrix Test") {

  // The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
  test_that("OptimatrixAdapter returns a optimatrix") {
    OptimatrixAdapterTestFixture fixture;
    bool result = fixture.TestOptimatrixReturnsCorrectValue(std::vector<int>(), std::vector<int>(), std::vector<double>());
    expect_true(result);
  }

}
