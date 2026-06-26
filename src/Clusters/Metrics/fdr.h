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

class FDR final : public ClusterMetric  {
    
public:
    FDR() : ClusterMetric("fdr") {};
    double getValue(const double tp, const double tn, const double fp, const double fn) override{
        const double pPrime = tp + fp;
        double falseDiscoveryRate = fp / pPrime;
        
        if(pPrime == 0)		{	falseDiscoveryRate = 0;		}
        
        if (std::isnan(falseDiscoveryRate) || std::isinf(falseDiscoveryRate)) { falseDiscoveryRate = 0; }
        
        return (1.0-falseDiscoveryRate);
    }
    std::string getCitation() override { return "http://www.mothur.org/wiki/FDR"; }
};

/***********************************************************************/



#endif /* fdr_hpp */
