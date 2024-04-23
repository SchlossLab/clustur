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

class TN : public ClusterMetric  {
    
public:
    TN() : ClusterMetric("tn") {};
    double getValue(double tp,  double tn,  double fp,  double fn){
                double tnmax = tn / (double)(tp + tn + fp + fn);
        
        if (std::isnan(tnmax) || std::isinf(tnmax)) { tnmax = 0; }
        
        return tnmax;
    }
    std::string getCitation() { return "http://www.mothur.org/wiki/TN"; }
    
private:
    
};

/***********************************************************************/




#endif /* tn_hpp */
