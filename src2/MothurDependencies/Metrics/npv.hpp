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

class NPV : public ClusterMetric  {
    
public:
    NPV() : ClusterMetric("npv") {};
    double getValue(double tp,  double tn,  double fp,  double fn){
        long long nPrime = tn + fn;
        double negativePredictiveValue = tn / (double) nPrime;
        
        if(nPrime == 0)		{	negativePredictiveValue = 0;		}
        
        if (std::isnan(negativePredictiveValue) || std::isinf(negativePredictiveValue)) { negativePredictiveValue = 0; }
        
        return negativePredictiveValue;
    }
    std::string getCitation() { return "http://www.mothur.org/wiki/NPV"; }
    
private:
    
};

/***********************************************************************/




#endif /* npv_hpp */
