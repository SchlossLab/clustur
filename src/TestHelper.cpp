//
// Created by Gregory Johnson on 4/15/24.
//

#include "TestHelpers/TestHelper.h"

void TestHelper::Print(const std::string &val) {
#if DEBUG_RCPP
    Rcpp::Rcout << val;
#else
    std::cout << val;
#endif

}
