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

class TP : public ClusterMetric  {
    
public:
    TP() : ClusterMetric("tp") {};
    double getValue(double tp,  double tn,  double fp,  double fn){
        double tpmax = tp / (double)(tp + tn + fp + fn);
        
        if (std::isnan(tpmax) || std::isinf(tpmax)) { tpmax = 0; }
        
        return tpmax;
    }
    std::string getCitation() { return "http://www.mothur.org/wiki/TP"; }
    
private:
    
};

/***********************************************************************/



#endif /* tp_hpp */
