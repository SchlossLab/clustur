//
//  sensspeccalc.hpp
//  Mothur
//
//  Created by Sarah Westcott on 1/22/18.
//  Copyright © 2018 Schloss Lab. All rights reserved.
//

#ifndef sensspeccalc_hpp
#define sensspeccalc_hpp

#include <vector>
#include "../DataStructures/ListVector.h"
#include "../DataStructures/OptiData.h"

class SensSpecCalc {

public:
    SensSpecCalc(const OptiData& matrix, const ListVector* list);
    ~SensSpecCalc() = default;

    void getResults(const OptiData& matrix, double& tp, double& tn, double& fp, double& fn) const;

private:
    std::vector<std::vector< int> > otus;
};



#endif /* sensspeccalc_hpp */
