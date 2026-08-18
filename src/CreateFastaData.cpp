//
// Created by Gregory Johnson on 8/18/26.
//
#include <Rcpp.h>
#include <string>
#include <vector>

#include "DataStructures/FastaDatabase.h"

//[[Rcpp::export]]
SEXP CreateFastaDataBase(const std::vector<std::string>& names, const std::vector<std::string>& sequences) {
    FastaDatabase* db = new FastaDatabase(names, sequences);
    return Rcpp::XPtr<FastaDatabase>(db);
}

//[[Rcpp::export]]
int GetSequences(const SEXP& fastaDatabase) {
    const Rcpp::XPtr<FastaDatabase> data(fastaDatabase);
    Rcpp::Rcout << data->GetFastaDataBase()[0].name << std::endl;
    Rcpp::Rcout << data->GetFastaDataBase()[0].sequence << std::endl;
    int res = data->GetFastaDataBase().size();
    return res;
}

