//
// Created by Gregory Johnson on 10/3/24.
//

#ifndef DISTANCEFILESTOSPARSE_H
#define DISTANCEFILESTOSPARSE_H
#include <string>
#include <vector>
#include "CountTableAdapter.h"
#include "../RowData.h"

class DistanceFilesToSparse {
public:
    DistanceFilesToSparse(CountTableAdapter  countTable, std::string method, std::string filePath);
    std::vector<std::vector<double>> CreateSparseList() const;
private:
    std::string filePath;
    std::string method;
    CountTableAdapter countTable;
};



#endif //DISTANCEFILESTOSPARSE_H
