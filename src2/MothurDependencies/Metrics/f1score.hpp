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

class F1Score : public ClusterMetric  {
    
public:
    F1Score() : ClusterMetric("f1score") {};
    double getValue(double tp,  double tn,  double fp,  double fn){
        long long p = 2.0 * tp;
        long long pPrime = fn + fp;
        double f1Score = 2.0 * tp / (double) (p + pPrime);
        
        if(p + pPrime == 0)	{	f1Score = 0;	}
        
        if (std::isnan(f1Score) || std::isinf(f1Score)) { f1Score = 0; }
        
        return f1Score;
    }
    std::string getCitation() { return "http://www.mothur.org/wiki/F1Score"; }
    
private:
    
};

/***********************************************************************/




#endif /* f1score_hpp */
