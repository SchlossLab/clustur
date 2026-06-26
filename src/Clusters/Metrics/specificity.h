//
//  specificity.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef specificity_hpp
#define specificity_hpp
#include <cmath>
#include "../ClusterMetric.h"
#include <string>


/***********************************************************************/

class Specificity final : public ClusterMetric  {
    
public:
    Specificity() : ClusterMetric("spec") {};
    double getValue(const double tp, const double tn, const double fp, const double fn) override {
        const double n = fp + tn;
        double specificity = tn / n;
        
        if(n == 0)			{	specificity = 0;	}
        if (std::isnan(specificity) || std::isinf(specificity)) { specificity = 0; }
        
        return specificity;
    } //ignores tp, fn
    std::string getCitation() override { return "http://www.mothur.org/wiki/Specificity"; }
};

/***********************************************************************/



#endif /* specificity_hpp */
