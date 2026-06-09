//
//  accuracy.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/11/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef accuracy_hpp
#define accuracy_hpp
#include <cmath>
#include "../ClusterMetric.h"
#include <string>


/***********************************************************************/

class Accuracy final : public ClusterMetric  {
    
public:
    Accuracy() : ClusterMetric("accuracy") {};
    double getValue(const double tp, const double tn, const double fp, const double fn) override
    {
        const double p = tp + fn;
        const double n = fp + tn;
        double accuracy = (tp + tn) / (p + n);
        if(p + n == 0) accuracy = 0;
        
        if (std::isnan(accuracy) || std::isinf(accuracy)) { accuracy = 0; }
        
        return accuracy;
    };
    std::string getCitation() override { return "http://www.mothur.org/wiki/Accuracy"; }
};

/***********************************************************************/


#endif /* accuracy_hpp */
