//
// Created by Gregory Johnson on 8/18/26.
//

#include "DataStructures/FastaDatabase.h"

FastaDatabase::FastaDatabase(const std::vector<std::string> &names,
    const std::vector<std::string> &sequences) {
    fastaDatabase.resize(names.size());
    for (size_t i = 0; i < names.size(); ++i) {
        fastaDatabase[i] = FastaData{names[i], sequences[i]};
    }
}
