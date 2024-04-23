//
//  accuracy.hpp
//  Mothur
//
//  Created by Sarah Westcott on 4/11/17.
//  Copyright © 2017 Schloss Lab. All rights reserved.
//

#ifndef accuracy_hpp
#define accuracy_hpp
#include <cmath>
#include "../ClusterMetric.h"
#include <string>


/***********************************************************************/

class Accuracy : public ClusterMetric  {
    
public:
    Accuracy() : ClusterMetric("accuracy") {};
    double getValue(double tp,  double tn,  double fp,  double fn)
    {
        long long p = tp + fn;
        long long n = fp + tn;
        double accuracy = (tp + tn) / (double) (p + n);
        if(p + n == 0)		{	accuracy = 0;								}
        
        if (std::isnan(accuracy) || std::isinf(accuracy)) { accuracy = 0; }
        
        return accuracy;
    };
    std::string getCitation() { return "http://www.mothur.org/wiki/Accuracy"; }
    
private:
    
};

/***********************************************************************/


#endif /* accuracy_hpp */
