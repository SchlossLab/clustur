//
// Created by Gregory Johnson on 8/21/26.
//

#include <Rcpp.h>
#include <fstream>

#include "Adapters/MatrixAdapter.h"
#include "FileReaders/ColumnDistanceMatrixReader.h"
#include "FileReaders/DistanceFileReader.h"
#include "FileReaders/ReadPhylipMatrix.h"

//[[Rcpp::export]]
bool DetermineIfPhylipOrColumnFile(const std::string& filePath) {
    std::ifstream data(filePath);

    if(!data.is_open()) {
        Rcpp::Rcout << "Please enter a valid file path\n";
        Rcpp::stop("Invalid file path.");
    }
    std::string line;
    std::getline(data, line);
    bool isPhylip = true;
    int count = 0;
    for(const auto& character: line) {
        if(character != '\n' && std::isspace(character))
            count++;
    }
    if(count > 1)
        isPhylip = false;
    data.close();
    return isPhylip;
}

//[[Rcpp::export]]
SEXP ProcessDistanceFiles(const std::string& filePath, const Rcpp::DataFrame& countTable, const double cutoff,
    const bool isSim) {
    const bool isPhylip = DetermineIfPhylipOrColumnFile(filePath);
    CountTableAdapter countTableAdapter;
    countTableAdapter.CreateDataFrameMap(countTable);
    if(isPhylip) {
        DistanceFileReader* read = new ReadPhylipMatrix(cutoff, isSim);
        read->SetCountTableAdapter(countTableAdapter);
        if (const bool succeeded = read->Read(filePath); !succeeded) {
            const std::set<std::string> params = read->GetFailureParameters();
            delete read;

            Utils::CheckForDistanceFileError(params);
        }

        return Rcpp::XPtr<DistanceFileReader>(read);
    }
    DistanceFileReader* read = new ColumnDistanceMatrixReader(cutoff, isSim);
    read->SetCountTableAdapter(countTableAdapter);
    if (const bool succeeded = read->Read(filePath); !succeeded) {
        const std::set<std::string> params = read->GetFailureParameters();
        delete read;
        Utils::CheckForDistanceFileError(params);
    }
    return Rcpp::XPtr<DistanceFileReader>(read);
}

 //[[Rcpp::export]]
 SEXP ProcessSparseMatrix(const std::vector<int> &xPosition,
     const std::vector<int> &yPosition, const std::vector<double> &data, const Rcpp::DataFrame& countTable,
     const double cutoff, const bool isSim) {
     CountTableAdapter countTableAdapter;
     countTableAdapter.CreateDataFrameMap(countTable);
     const MatrixAdapter adapter(xPosition, yPosition, data, cutoff, isSim, countTableAdapter);
     auto* read = new DistanceFileReader(adapter.CreateSparseMatrix(),
         adapter.CreateListVector(), countTableAdapter, cutoff, isSim);
     return Rcpp::XPtr<DistanceFileReader>(read);
 }


//[[Rcpp::export]]
double AddDataToDistanceData(SEXP& refData, const SEXP& fitData,
    const Rcpp::DataFrame& refFasta,
    const Rcpp::DataFrame& fitFasta,
    const double cutoff) {

    const Rcpp::XPtr<DistanceFileReader> refDistanceData(refData);
    const CountTableAdapter refCountTableAdapter = refDistanceData.get()->GetCountTableAdapter();
    const FastaDatabase refFastaDatabase(refFasta["sequence_name"], refFasta["sequence"]);

    const Rcpp::XPtr<DistanceFileReader> fitDistanceData(fitData);
    const CountTableAdapter fitCountTableAdapter = fitDistanceData.get()->GetCountTableAdapter();
    const SparseDistanceMatrix* fitSparseMatrix =  fitDistanceData.get()->GetSparseMatrix();
    const ListVector* fitListVector = fitDistanceData.get()->GetListVector();
    FastaDatabase fitFastaDatabase(fitFasta["sequence_name"], fitFasta["sequence"]);

    refDistanceData.get()->AddFittedDataToReference(fitSparseMatrix, fitListVector,
        fitCountTableAdapter, refFastaDatabase, fitFastaDatabase, cutoff);
    const double fitCountSize = static_cast<double>(fitCountTableAdapter.GetSequences().size());
    const double refCountSize = static_cast<double>(refCountTableAdapter.GetSequences().size());
    const double result = fitCountSize / (fitCountSize + refCountSize);
    return result;// Fit Percentage
}

//[[Rcpp::export]]
SEXP CopyObject(const SEXP& distanceObject) {
    const Rcpp::XPtr<DistanceFileReader> refDistanceData(distanceObject);
    DistanceFileReader* copy = new DistanceFileReader(*refDistanceData.get());
    return Rcpp::XPtr<DistanceFileReader>(copy);
}
