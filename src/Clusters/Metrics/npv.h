//
//  npv.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/11/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef npv_hpp
#define npv_hpp
#include <cmath>
#include "../ClusterMetric.h"
#include <string>


/***********************************************************************/

class NPV final : public ClusterMetric  {
    
public:
    NPV() : ClusterMetric("npv") {};
    double getValue(const double tp, const double tn, const double fp, const double fn) override{
        const double nPrime = tn + fn;
        double negativePredictiveValue = tn / nPrime;
        
        if(nPrime == 0)	negativePredictiveValue = 0;
        
        if (std::isnan(negativePredictiveValue) || std::isinf(negativePredictiveValue)) { negativePredictiveValue = 0; }
        
        return negativePredictiveValue;
    }
    std::string getCitation() override { return "http://www.mothur.org/wiki/NPV"; }
};

/***********************************************************************/




#endif /* npv_hpp */
