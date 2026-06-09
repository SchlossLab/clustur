//
//  mcc.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef mcc_hpp
#define mcc_hpp
#include <cmath>
#include "../ClusterMetric.h"
#include <string>




/***********************************************************************/

class MCC final : public ClusterMetric  {
    
public:
    MCC() : ClusterMetric("mcc") {};
    double getValue(const double tp, const double tn, const double fp, const double fn) override {
        const double p = tp + fn;
        const double n = fp + tn;
        const double pPrime = tp + fp;
        const double nPrime = tn + fn;
        
        double matthewsCorrCoef = ((tp * tn) - (fp * fn)) / sqrt(p * n * pPrime * nPrime);
        
        if(p == 0 || n == 0 || pPrime == 0 || nPrime == 0){	matthewsCorrCoef = 0;	}
        
        if (std::isnan(matthewsCorrCoef) || std::isinf(matthewsCorrCoef)) { matthewsCorrCoef = 0; }
        
        return matthewsCorrCoef;
    }
    std::string getCitation() override { return "http://www.mothur.org/wiki/MCC"; }
};

/***********************************************************************/

#endif /* mcc_hpp */
