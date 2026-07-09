//
// Created by Gregory Johnson on 7/9/26.
//

#ifndef REFACTOR_SUBSAMPLE_H
#define REFACTOR_SUBSAMPLE_H


//
//  subsample.h
//  Mothur
//
//  Created by Sarah Westcott on 4/2/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include <map>
#include <set>

//subsampling overwrites the sharedRabunds.  If you need to reuse the original use the getSamplePreserve function.

class SubSample {

    public:
        SubSample() = default;
        ~SubSample() = default;
        static std::set<long long> getWeightedSample(std::map<long long, long long>&, long long); //map of sequence names -> weight (could be abundance or some other measure), num to sample

};

#endif //REFACTOR_SUBSAMPLE_H