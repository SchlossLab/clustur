//
// Created by Gregory Johnson on 9/2/26.
//

#include <utility>

#include "Clusters/ClusterSplit.h"

#include "Clusters/OptiCluster.h"
#include "Adapters/OptimatrixAdapter.h"
#include "DataStructures/OptiData.h"
#include "DataStructures/SplitMatrix.h"

ClusterSplit::ClusterSplit(FastaDatabase fastaDatabase, const std::vector<TaxonomyData> &taxaData,
	PairwiseDistanceCalculator *calculator,
	ClusterParameters parameters,
	ClusterMetric *metric,
	CountTableAdapter adapter,
	const double cutoff,
	const int taxonomyCutoff,
	const int seed,
	const int numberOfThreads):fastaData(std::move(fastaDatabase)), calculator(calculator),
		clusterParameters(parameters), metric(metric), taxaData(taxaData), countTableAdapter(std::move(adapter)),
		taxonomyCutoff(taxonomyCutoff), seed(seed), numberOfThreads(numberOfThreads), cutoff(cutoff){}

ClusterExport * ClusterSplit::Execute() {
	std::vector<OptiDataComponent> splitMatrices =
		SplitMatrix::splitClassify(taxaData, fastaData, calculator, cutoff, taxonomyCutoff, numberOfThreads);

	//****************** break up files between processes and cluster each file set ******************************//
	std::set<std::string> labels;
    const std::vector<ClusterExport*> results = createProcesses(splitMatrices, labels, numberOfThreads);

	ClusterExport* completeList = mergeLists(results);
	// // Add singletons
	const std::vector<std::string>& allSequences = countTableAdapter.GetSequences();
	const ListVector& completedListVector = completeList->GetListVector().listVector;
	ListVector singletons(static_cast<int>(allSequences.size() - completedListVector.getNumSeqs()));
	std::unordered_set<std::string> nonSingletons;

	nonSingletons.reserve(completedListVector.size());

	for (int i = 0; i < completedListVector.size(); i++) {
		std::vector<std::string> splitString;
		Utils::splitAtComma(completedListVector.get(i), splitString);
		for (const auto& str : splitString) {
			nonSingletons.insert(str);
		}
	}

	int counter = 0;
	for (const auto& seq : allSequences) {
		if (nonSingletons.find(seq) != nonSingletons.end()) continue; // its a singleton
		singletons.set(counter++, seq);
	}
	// Add all the other sequences to singleton
	const int singletonSize = singletons.size();
	const int completedListVectorSize = completedListVector.size();
	singletons.resize(singletonSize + completedListVectorSize);

	for (int i = singletonSize; i < completedListVectorSize + singletonSize; i++) {
		singletons.set(i, completedListVector.get(i - singletonSize));
	}

	// singletons.push_back(completedListVector);
	completeList->SetListVector(singletons, std::to_string(cutoff));
	return completeList;
}
//**********************************************************************************************************************
ClusterExport* ClusterSplit::mergeLists(const std::vector<ClusterExport*>& exportedResults){

	size_t completeSize = 0;
	for (const auto& result : exportedResults) {
		completeSize += result->GetListVector().listVector.size();
	}
	ListVector completeListVector(completeSize);
	int count = 0;
	for (const auto& result : exportedResults) {
		const ListVector& listVector = result->GetListVector().listVector;
		for (long long i = 0; i < static_cast<long long>(listVector.size()); i++) {
			completeListVector.set(count++, listVector.get(i));
		}
	}
	ClusterExport* result = new ClusterExport();
	result->SetListVector(completeListVector, std::to_string(cutoff));
	return result;
}


std::vector<ClusterExport *> ClusterSplit::createProcesses(std::vector<OptiDataComponent> &distanceMatrices,
                                                           std::set<std::string> &labels,
                                                           size_t processors){
    //sanity check
	const size_t matricesSize = distanceMatrices.size();
    if (processors > matricesSize) {
	    processors = matricesSize;
    }
  //  deleteFiles = false; //so if we need to recalc the processors the files are still there
   // vector<string> listFiles;
    std::vector<std::vector<OptiDataComponent>> dividedWork(processors); //distNames[1] = vector of filenames for process 1...
    // dividedNames.resize(processors);

    //for each file group figure out which process will complete it
    //want to divide the load intelligently so the big files are spread between processes
    for (int i = 0; i < matricesSize; i++) {
        int processToAssign = (i+1) % processors;
        if (processToAssign == 0) {
	        processToAssign = processors;
        }

        dividedWork[(processToAssign-1)].emplace_back(distanceMatrices[i]);
        // if ((processToAssign-1) == 1) { m->mothurOut(distName[i].begin()->first + "\n"); }
    }

    //now lets reverse the order of ever other process, so we balance big files running with little ones
    for (int i = 0; i < processors; i++) {
        int remainder = ((i+1) % processors);
        if (remainder) { std::reverse(dividedWork[i].begin(), dividedWork[i].end());  }
    }

    // if (m->getControl_pressed()) { return listFiles; }


    //create array of worker threads
    std::vector<RcppThread::Thread*> workerThreads;
    std::vector<SplitClusterData*> data;
	data.reserve(dividedWork.size() - 1);

    //Lauch worker threads
    for (int i = 1; i < processors; i++) {
    	SplitClusterData* dataBundle = new SplitClusterData(dividedWork[i], clusterParameters, metric,
    		RandomNumberSitmo(seed + i), cutoff, i);
        data.emplace_back(dataBundle);

        workerThreads.push_back(new RcppThread::Thread([&, dataBundle] {
	        cluster(dataBundle);
        }));
    }


	SplitClusterData* dataBundle = new SplitClusterData(dividedWork[0], clusterParameters, metric,
		RandomNumberSitmo(seed), cutoff, 0);
    cluster(dataBundle);

	std::vector<ClusterExport*> results = dataBundle->results;
	results.reserve(distanceMatrices.size());
    for (int i = 0; i < processors - 1; i++) {
        workerThreads[i]->join();
        delete workerThreads[i];
    }
	for (const auto& result: data) {
		results.insert(results.end(), result->results.cbegin(), result->results.cend());
	}

	dataBundle->results.clear();
	delete dataBundle;
	for (auto & dat : data) {
		delete dat;
	}

    //deleteFiles = true;

    return results;
}

void ClusterSplit::cluster(SplitClusterData* params) const {
	params->results.reserve(params->dividedData.size());
	for (auto&[matrix, listVector] : params->dividedData) {
		if (params->clusterParameters.GetClusterType() == "opti") {
			OptimatrixAdapter adapter(cutoff);
			OptiData* data = adapter.ConvertToOptimatrix(&matrix, &listVector, false);
			ClusterMethod* method = new OptiCluster(data, params->metric, params->rng, cutoff, 0);
			method->SetClusterParameters(params->clusterParameters);
			params->results.emplace_back(method->Execute());
			delete method;
			delete data;
			continue;
		}
		RAbundVector rAbund = listVector.getRAbundVector();
		ClusterMethod* method = Utils::GetClusterMethod(params->clusterParameters.GetClusterType(), &listVector, &matrix,
			rAbund, cutoff);
		method->SetClusterParameters(params->clusterParameters);
		params->results.emplace_back(method->Execute());
		delete method;

	}
}