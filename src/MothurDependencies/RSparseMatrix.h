//
// Created by Gregory Johnson on 7/31/24.
//

#ifndef RSPARSEMATRIX_H
#define RSPARSEMATRIX_H
#include <vector>


struct RSparseData {
    RSparseData(const int i, const int j, const double value) {
        iIndex = i;
        jIndex = j;
        data = value;
    }
    int iIndex;
    int jIndex;
    double data;
};

struct RSpraseMatrix {
    RSpraseMatrix(const std::vector<int> &iIndex, const std::vector<int> &jIndex, const std::vector<double> &data){
        for(size_t i = 0; i < data.size(); i++) {
            sparseValues.emplace_back(iIndex[i], jIndex[i], data[i]);
        }
    }
    std::vector<RSparseData> sparseValues;
};


#endif //RSPARSEMATRIX_H