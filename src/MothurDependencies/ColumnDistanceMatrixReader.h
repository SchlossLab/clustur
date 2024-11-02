//
// Created by Gregory Johnson on 9/13/24.
//

#ifndef COLUMNDISTANCEMATRIXREADER_H
#define COLUMNDISTANCEMATRIXREADER_H
#include <string>
#include <fstream>
#include <map>
#include <vector>
#include "ListVector.h"
#include "SparseDistanceMatrix.h"
#include "../RowData.h"
#include "../Adapters/CountTableAdapter.h"
#include "../Adapters/DistanceFileReader.h"

class ColumnDistanceMatrixReader final : public DistanceFileReader {
public:
    ColumnDistanceMatrixReader(double cutoff, bool isSimularity);
    bool Read(const std::string &filePath) override;
};



#endif //COLUMNDISTANCEMATRIXREADER_H
