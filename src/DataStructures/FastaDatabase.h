//
// Created by Gregory Johnson on 8/18/26.
//

#ifndef REFACTOR_FASTA_H
#define REFACTOR_FASTA_H
#include <string>
#include <vector>


struct FastaData {
    std::string name;
    std::string sequence;
};

class FastaDatabase {
public:
    FastaDatabase(const std::vector<std::string>& names, const std::vector<std::string>& sequences);
    const std::vector<FastaData>& GetFastaDataBase() const {return fastaDatabase;}
private:
    std::vector<FastaData> fastaDatabase;
};


#endif //REFACTOR_FASTA_H