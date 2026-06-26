//
//  fn.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef fn_hpp
#define fn_hpp
#include <cmath>
#include "../ClusterMetric.h"
#include <string>
/***********************************************************************/

class FN : public ClusterMetric  {
    
public:
    FN() : ClusterMetric("fn") {};
    double getValue(double tp,  double tn,  double fp,  double fn){
        double fnmin = fn / (double)(tp + tn + fp + fn);
        
        if (std::isnan(fnmin) || std::isinf(fnmin)) { fnmin = 0; }
        
        return (1.0 - fnmin);
    }
    std::string getCitation() { return "http://www.mothur.org/wiki/FN"; }
    
private:
    
};

/***********************************************************************/



#endif /* fn_hpp */
