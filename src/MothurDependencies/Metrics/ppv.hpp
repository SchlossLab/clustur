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

class PPV : public ClusterMetric  {
    
public:
    PPV() : ClusterMetric("ppv") {};
    double getValue(double tp,  double tn,  double fp,  double fn){
        long long pPrime = tp + fp;
        double positivePredictiveValue = tp / (double) pPrime;
        
        if(pPrime == 0)		{	positivePredictiveValue = 0;		}
        
        if (std::isnan(positivePredictiveValue) || std::isinf(positivePredictiveValue)) { positivePredictiveValue = 0; }
        return positivePredictiveValue;
        
    }
    std::string getCitation() { return "http://www.mothur.org/wiki/PPV"; }
    
private:
    
};

/***********************************************************************/



#endif /* ppv_hpp */
