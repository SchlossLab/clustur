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

class FPFN final : public ClusterMetric  {
    
public:
    FPFN() : ClusterMetric("fpfn") {};
    double getValue(const double tp, const double tn, const double fp, const double fn) override{
        const double p = fp + fn;
        
        double fpfn = 1.0 - (p / (tp + tn + fp + fn)); //minimize
        
        if (std::isnan(fpfn) || std::isinf(fpfn)) { fpfn = 0; }
        
        return fpfn;
    }
    std::string getCitation() override { return "http://www.mothur.org/wiki/FPFN"; }
};

/***********************************************************************/


#endif /* fpfn_hpp */
