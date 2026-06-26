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

class Sensitivity : public ClusterMetric  {
    
public:
    Sensitivity() : ClusterMetric("sens") {};
    double getValue(double tp,  double tn,  double fp,  double fn){
        long long p = tp + fn;
        double sensitivity = tp / (double) p;
        
        if(p == 0)	{	sensitivity = 0;	}
        if (std::isnan(sensitivity) || std::isinf(sensitivity)) { sensitivity = 0; }
        
        return sensitivity;
    } //ignores tn, fp
    std::string getCitation() { return "http://www.mothur.org/wiki/Sensitivity"; }
    
private:
    
};

/***********************************************************************/

#endif /* sensitivity_hpp */
