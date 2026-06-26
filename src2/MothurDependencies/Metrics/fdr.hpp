//
//  fdr.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/11/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef fdr_hpp
#define fdr_hpp
#include <cmath>
#include "../ClusterMetric.h"
#include <string>

/***********************************************************************/

class FDR : public ClusterMetric  {
    
public:
    FDR() : ClusterMetric("fdr") {};
    double getValue(double tp,  double tn,  double fp,  double fn){
        long long pPrime = tp + fp;
        double falseDiscoveryRate = fp / (double) pPrime;
        
        if(pPrime == 0)		{	falseDiscoveryRate = 0;		}
        
        if (std::isnan(falseDiscoveryRate) || std::isinf(falseDiscoveryRate)) { falseDiscoveryRate = 0; }
        
        return (1.0-falseDiscoveryRate);
    }
    std::string getCitation() { return "http://www.mothur.org/wiki/FDR"; }
    
private:
    
};

/***********************************************************************/



#endif /* fdr_hpp */
