//
//  sensspeccalc.cpp
//  Mothur
//
//  Created by Sarah Westcott on 1/22/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#include "MothurDependencies/SensSpecCalc.h"
#include "MothurDependencies/Utils.h"
#include <string>
#include <map>


//***************************************************************************************************************
//removes anyone with no valid dists and changes name to matrix short names
SensSpecCalc::SensSpecCalc(const OptiData& matrix, const ListVector* list){
        std::map<std::string, long long> nameIndex = matrix.getNameIndexMap();

        if (list != nullptr) {
            //for each bin
            for (int i = 0; i < list->getNumBins(); i++) {

                std::string binnames = list->get(i);
                std::vector<std::string> bnames;
                Utils::splitAtComma(binnames, bnames);

                std::vector<int> newNames;
                for (const auto & name : bnames) {
                    auto itSeq1 = nameIndex.find(name);
                    long long seq1Index = -1;
                    if (itSeq1 != nameIndex.end()) { seq1Index = itSeq1->second; } //you have distances in the matrix

                    newNames.push_back(seq1Index);
                }

                //if there are names in this bin add to new list
                if (!newNames.empty()) { otus.push_back(newNames); }
            }
        }
}
//***************************************************************************************************************
void SensSpecCalc::getResults(const OptiData& matrix, double& tp, double& tn, double& fp, double& fn) const {
        tp = 0; tn = 0; fp = 0; fn = 0;

        for(auto & otu : otus){

            for(int i=0;i<otu.size();i++){
                for(int j=0;j<i;j++){
                    if (matrix.isClose(otu[i], otu[j])) { tp++; }
                    else { fp++; }
                }
            }
        }
        const double numSeqs = static_cast<double>(matrix.getNumSeqs() + matrix.getNumSingletons());
        const double numDists = matrix.getNumDists(); //square matrix OptiData:: uses the parent class function so that we can pass a optiref matrix

        fn = (numDists/2) - tp;
        tn = numSeqs * (numSeqs-1)/2  - (fp + fn + tp);
}

//***************************************************************************************************************

