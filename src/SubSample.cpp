//
// Created by Gregory Johnson on 7/9/26.
//
#include <Rcpp.h>
#include "MothurDependencies/Utils.h"
#include "RNG/SubSample.h"

std::set<long long> SubSample::getWeightedSample(std::map<long long, long long>& nameMap, const long long num) {
    std::set<long long> sampleNames;

    const long long totalSeqs = nameMap.size();
    if (totalSeqs < num) { Rcpp::stop("[ERROR]: Requesting sample size larger than number of seqeunces, quitting.\n"); }
    if (totalSeqs == num) {
        for (auto &[fst, snd] : nameMap) { sampleNames.insert(fst);     }
        return sampleNames;
    }
    
    long long numSampled = 0;
    std::map<long long, std::set<long long> > weights;//weight -> names of seqs with that weight
    std::map<long long, std::set<long long> >::iterator itWeight;
    long long total = 0;
    for (auto &[first, second] : nameMap) {
        total += second;
        itWeight = weights.find(second);
        if (itWeight == weights.end()) { //this is a weight we haven't seen before
            std::set<long long> temp;
            temp.insert(first);
            weights[second] = temp;
        }else {
            weights[second].insert(first); //dup weight, combine to save memory
        }
    }
    
    //find running total
    long long runningTotal = 0;
    std::map<long long, std::set<long long> > cumulative;//weight + sum so far -> names of seqs with that weight
    for (itWeight = weights.begin(); itWeight != weights.end(); itWeight++) {
        const size_t count = itWeight->second.size(); //number of seqs with this weight
        runningTotal += itWeight->first * count;
        cumulative[runningTotal] = itWeight->second;
    }
    weights.clear();
    
    while(numSampled != num) {
        long long index = Utils::getRandomIndex(total); //random index including weights
        
        std::map<long long, std::set<long long> >::iterator itWeight = cumulative.lower_bound(index);
        
        sampleNames.insert(*itWeight->second.begin()); //save name in sample names
        
        itWeight->second.erase(itWeight->second.begin()); //remove seq since we sampled it
        
        if (itWeight->second.empty()) { cumulative.erase(itWeight); total = cumulative.rbegin()->first;  } //remove this weight if all seqs are sampled. Reset bound if needed.
        
        numSampled = sampleNames.size();
    }

    return sampleNames;
}
