//
//  tn.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef tn_hpp
#define tn_hpp
#include <cmath>
#include "../ClusterMetric.h"
#include <string>


/***********************************************************************/

class TN final : public ClusterMetric  {
    
public:
    TN() : ClusterMetric("tn") {};
    double getValue(const double tp, const double tn, const double fp, const double fn) override{
        double tnmax = tn / (tp + tn + fp + fn);
        
        if (std::isnan(tnmax) || std::isinf(tnmax)) { tnmax = 0; }
        
        return tnmax;
    }
    std::string getCitation() override { return "http://www.mothur.org/wiki/TN"; }
};

/***********************************************************************/




#endif /* tn_hpp */
