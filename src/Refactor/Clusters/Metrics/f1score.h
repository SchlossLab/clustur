//
//  f1score.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/11/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef f1score_hpp
#define f1score_hpp
#include <cmath>
#include "../ClusterMetric.h"
#include <string>


/***********************************************************************/

class F1Score final : public ClusterMetric  {
    
public:
    F1Score() : ClusterMetric("f1score") {};
    double getValue(const double tp, const double tn, const double fp, const double fn) override{
        const double p = 2.0 * tp;
        const double pPrime = fn + fp;
        double f1Score = 2.0 * tp / (p + pPrime);
        
        if(p + pPrime == 0)	{	f1Score = 0;	}
        
        if (std::isnan(f1Score) || std::isinf(f1Score)) { f1Score = 0; }
        
        return f1Score;
    }
    std::string getCitation() override { return "http://www.mothur.org/wiki/F1Score"; }
};

/***********************************************************************/




#endif /* f1score_hpp */
