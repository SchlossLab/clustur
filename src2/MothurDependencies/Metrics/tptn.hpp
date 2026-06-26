//
//  tptn.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/10/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef tptn_hpp
#define tptn_hpp
#include <cmath>
#include "../ClusterMetric.h"
#include <string>


/***********************************************************************/

class TPTN : public ClusterMetric  {
    
public:
    TPTN() : ClusterMetric("tptn") {};
        double getValue(double tp,  double tn,  double fp,  double fn)  {
        long long p = tp + tn;
        double tptn = p / (double)(tp + tn + fp + fn);
        
        if (std::isnan(tptn) || std::isinf(tptn)) { tptn = 0; }
        return tptn;

    }   
    std::string getCitation() { return "http://www.mothur.org/wiki/TPTN"; }
    
private:
    
};

/***********************************************************************/




#endif /* tptn_hpp */
