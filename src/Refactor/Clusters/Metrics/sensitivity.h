//
//  sensitivity.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef sensitivity_hpp
#define sensitivity_hpp
#include <cmath>
#include "../ClusterMetric.h"
#include <string>


/***********************************************************************/

class Sensitivity final : public ClusterMetric  {
    
public:
    Sensitivity() : ClusterMetric("sens") {};
    double getValue(const double tp, const double tn, const double fp, const double fn) override {
        const double p = tp + fn;
        double sensitivity = tp / p;
        
        if(p == 0)	{sensitivity = 0; }
        if (std::isnan(sensitivity) || std::isinf(sensitivity)) { sensitivity = 0; }
        
        return sensitivity;
    } //ignores tn, fp
    std::string getCitation() override { return "http://www.mothur.org/wiki/Sensitivity"; }
};

/***********************************************************************/

#endif /* sensitivity_hpp */
