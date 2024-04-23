//
//  fpfn.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef fpfn_hpp
#define fpfn_hpp
#include <cmath>
#include "../ClusterMetric.h"
#include <string>

/***********************************************************************/

class FPFN : public ClusterMetric  {
    
public:
    FPFN() : ClusterMetric("fpfn") {};
    double getValue(double tp,  double tn,  double fp,  double fn){
        long long p = fp + fn;
        
        double fpfn = 1.0 - (p / (double)(tp + tn + fp + fn)); //minimize
        
        if (std::isnan(fpfn) || std::isinf(fpfn)) { fpfn = 0; }
        
        return fpfn;
    }
    std::string getCitation() { return "http://www.mothur.org/wiki/FPFN"; }
    
private:
    
};

/***********************************************************************/


#endif /* fpfn_hpp */
