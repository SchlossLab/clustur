#include <Rcpp.h>
#include <cctype>
#include <algorithm>

#include "Adapters/CountTableAdapter.h"
#include "Adapters/OptimatrixAdapter.h"
#include "SharedFileData/SharedFileBuilder.h"
#include "FileReaders/DistanceFileReader.h"
#include "MothurDependencies/CreateDistanceMatrix.h"
#include "MothurDependencies/OneGapPairwiseDistance.h"
#include "MothurDependencies/PairwiseDistanceCalculator.h"


//[[Rcpp::export]]
Rcpp::DataFrame GetDistanceDataFrame(const SEXP& fileReader) {
    const Rcpp::XPtr<DistanceFileReader> ptr(fileReader);
    return ptr.get()->SparseMatrixToDataFrame();
}

//[[Rcpp::export]]
Rcpp::DataFrame GetCountTable(const SEXP& fileReader) {
    const Rcpp::XPtr<DistanceFileReader> ptr(fileReader);
    return ptr.get()->GetCountTable();
}


//[[Rcpp::export]]
Rcpp::DataFrame CreateDataFrameFromSparseCountTable(const Rcpp::DataFrame& countTable) {
    CountTableAdapter adapter;
    adapter.CreateDataFrameMapFromSparseCountTable(countTable);
    return adapter.ReCreateDataFrame();
}


//[[Rcpp::export]]
double CreateSparseMatrix(const std::vector<std::string>& sequences, const double cutoff) {
    PairwiseDistanceCalculator* calculator = new OneGapPairwiseDistance();
    SparseDistanceMatrix matrix = CreateDistanceMatrix::CreateSparseDistanceMatrix(sequences, calculator, cutoff);
    delete calculator;
    return 0;

}
