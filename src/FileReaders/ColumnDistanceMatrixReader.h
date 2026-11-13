//
// Created by Gregory Johnson on 9/13/24.
//

#ifndef COLUMNDISTANCEMATRIXREADER_H
#define COLUMNDISTANCEMATRIXREADER_H
#include "DistanceFileReader.h"

class ColumnDistanceMatrixReader final : public DistanceFileReader {
public:
    ColumnDistanceMatrixReader(double cutoff, bool isSimularity);
    ~ColumnDistanceMatrixReader() override = default;
    bool Read(const std::string &filePath) override;
};



#endif //COLUMNDISTANCEMATRIXREADER_H
