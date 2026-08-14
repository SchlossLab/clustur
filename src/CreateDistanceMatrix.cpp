//
// Created by Gregory Johnson on 8/13/26.
//

#include "MothurDependencies/CreateDistanceMatrix.h"
#include <algorithm>

SparseDistanceMatrix CreateDistanceMatrix::CreateSparseDistanceMatrix(const std::vector<std::string> &sequences,
                                                                      PairwiseDistanceCalculator* calculator,
                                                                      const double cutoff) {
    const size_t sequenceLength = sequences[0].size();
    const bool res = std::all_of(sequences.begin(), sequences.end(), [sequenceLength](const std::string &sequence) {
        return sequenceLength == sequence.size();
    });
    if (!res) Rcpp::stop("Not all sequences are of equal length");

    const size_t sequenceSize = sequences.size();

    SparseDistanceMatrix sparseMatrix(sequenceSize);
    for (size_t i = 0; i < sequenceSize - 1; i++) {
        for (size_t j = i + 1; j < sequenceSize; j++) {
            if (const double dist = calculator->Execute(sequences[i], sequences[j]); dist < cutoff) {
                const PDistCell value(i, dist);
                sparseMatrix.addCell(j, value);
            }
        }
    }

    return sparseMatrix;
}
