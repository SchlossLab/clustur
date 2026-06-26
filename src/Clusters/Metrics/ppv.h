//
//  ppv.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/11/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef ppv_hpp
#define ppv_hpp
#include <cmath>
#include "../ClusterMetric.h"
#include <string>


/***********************************************************************/

class PPV final : public ClusterMetric  {
    
public:
    PPV() : ClusterMetric("ppv") {};
    double getValue(const double tp, const double tn, const double fp, const double fn) override{
        const double pPrime = tp + fp;
        double positivePredictiveValue = tp / pPrime;
        
        if(pPrime == 0)	positivePredictiveValue = 0;
        
        if (std::isnan(positivePredictiveValue) || std::isinf(positivePredictiveValue)) { positivePredictiveValue = 0; }
        return positivePredictiveValue;
        
    }
    std::string getCitation() override { return "http://www.mothur.org/wiki/PPV"; }
};

/***********************************************************************/



#endif /* ppv_hpp */
