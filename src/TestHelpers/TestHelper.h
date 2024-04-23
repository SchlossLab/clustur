//
// Created by Gregory Johnson on 4/15/24.
//

#ifndef TESTHELPER_H
#define TESTHELPER_H
#define DEBUG_RCPP 1
#if DEBUG_RCPP
    #include <Rcpp.h>
#endif
#include <string>
#include <iostream>

class TestHelper {
public:
    static void Print(const std::string&);
};



#endif //TESTHELPER_H
