//
//  tp.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef tp_hpp
#define tp_hpp
#include <cmath>
#include "../ClusterMetric.h"
#include <string>


/***********************************************************************/

class TP final : public ClusterMetric  {
    
public:
    TP() : ClusterMetric("tp") {};
    double getValue(const double tp, const double tn, const double fp, const double fn) override{
        double tpmax = tp / (tp + tn + fp + fn);
        
        if (std::isnan(tpmax) || std::isinf(tpmax)) { tpmax = 0; }
        
        return tpmax;
    }
    std::string getCitation() override { return "http://www.mothur.org/wiki/TP"; }
};

/***********************************************************************/



#endif /* tp_hpp */
