//
//  fp.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef fp_hpp
#define fp_hpp
#include <cmath>
#include "../ClusterMetric.h"
#include <string>


/***********************************************************************/

class FP : public ClusterMetric  {
    
public:
    FP() : ClusterMetric("fp") {};
    double getValue(double tp,  double tn,  double fp,  double fn){
        double fpmin = fp / (double)(tp + tn + fp + fn);
        
        if (std::isnan(fpmin) || std::isinf(fpmin)) { fpmin = 0; }
        
        return (1.0 - fpmin);
    }
    std::string getCitation() { return "http://www.mothur.org/wiki/FP"; }
    
private:
    
};

/***********************************************************************/



#endif /* fp_hpp */
