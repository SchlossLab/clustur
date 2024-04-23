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

class MCC : public ClusterMetric  {
    
public:
    MCC() : ClusterMetric("mcc") {};
    double getValue(double tp,  double tn,  double fp,  double fn){
        double p = tp + fn;
        double n = fp + tn;
        double pPrime = tp + fp;
        double nPrime = tn + fn;
        
        double matthewsCorrCoef = ((tp * tn) - (fp * fn)) / sqrt(p * n * pPrime * nPrime);
        
        if(p == 0 || n == 0 || pPrime == 0 || nPrime == 0){	matthewsCorrCoef = 0;	}
        
        if (std::isnan(matthewsCorrCoef) || std::isinf(matthewsCorrCoef)) { matthewsCorrCoef = 0; }
        
        return matthewsCorrCoef;
    }
    std::string getCitation() { return "http://www.mothur.org/wiki/MCC"; }
    
private:
    
};

/***********************************************************************/

#endif /* mcc_hpp */
