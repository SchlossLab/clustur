//
// Created by Gregory Johnson on 9/2/26.
//

#include <utility>

#include "Clusters/ClusterSplit.h"

#include "Clusters/OptiCluster.h"
#include "Adapters/OptimatrixAdapter.h"
#include "Clusters/Metrics/accuracy.h"
#include "Clusters/Metrics/f1score.h"
#include "Clusters/Metrics/fdr.h"
#include "Clusters/Metrics/mcc.h"
#include "Clusters/Metrics/npv.h"
#include "Clusters/Metrics/ppv.h"
#include "Clusters/Metrics/sensitivity.h"
#include "Clusters/Metrics/specificity.h"
#include "DataStructures/OptiData.h"
#include "DataStructures/SplitMatrix.h"
#include "MothurDependencies/SensSpecCalc.h"

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
	if (singletons.size() == 0) {
		return completeList;
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
	if (clusterParameters.GetClusterType() == "opti") RunSenSpec();
	//delete results
	for (const auto& result : results) {
		delete result;
	}
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
    for (size_t i = 0; i < matricesSize; i++) {
        size_t processToAssign = (i+1) % processors;
        if (processToAssign == 0) {
	        processToAssign = processors;
        }

        dividedWork[(processToAssign-1)].emplace_back(distanceMatrices[i]);
        // if ((processToAssign-1) == 1) { m->mothurOut(distName[i].begin()->first + "\n"); }
    }

    //now lets reverse the order of ever other process, so we balance big files running with little ones
    for (size_t i = 0; i < processors; i++) {
        const size_t remainder = ((i+1) % processors);
        if (remainder) { std::reverse(dividedWork[i].begin(), dividedWork[i].end());  }
    }

    // if (m->getControl_pressed()) { return listFiles; }


    //create array of worker threads
    std::vector<RcppThread::Thread*> workerThreads;
    std::vector<SplitClusterData*> data;
	data.reserve(dividedWork.size() - 1);

    //Lauch worker threads
    for (int i = 1; i < static_cast<int>(processors); i++) {
    	SplitClusterData* dataBundle = new SplitClusterData(dividedWork[i], clusterParameters, metric,
    		seed + i, cutoff);
        data.emplace_back(dataBundle);

        workerThreads.push_back(new RcppThread::Thread([&, dataBundle] {
	        cluster(dataBundle);
        }));
    }


	SplitClusterData* dataBundle = new SplitClusterData(dividedWork[0], clusterParameters, metric,
		seed, cutoff);
    cluster(dataBundle);

	std::vector<ClusterExport*> results = dataBundle->results;
	results.reserve(distanceMatrices.size());
	if (dataBundle->clusterParameters.GetClusterType() == "opti") {
		fp += dataBundle->fp;
		fn += dataBundle->fn;
		tp += dataBundle->tp;
		tn += dataBundle->tn;
	}
    for (size_t i = 0; i < processors - 1; i++) {
        workerThreads[i]->join();
        delete workerThreads[i];
    }
	for (const auto& result: data) {
		results.insert(results.end(), result->results.cbegin(), result->results.cend());
		if (result->clusterParameters.GetClusterType() != "opti") continue;
		fp += result->fp;
		fn += result->fn;
		tp += result->tp;
		tn += result->tn;
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
			double fp = 0, fn = 0, tp = 0, tn = 0;
			OptimatrixAdapter adapter(cutoff);
			OptiData* data = adapter.ConvertToOptimatrix(&matrix, &listVector, false);
			ClusterMethod* method = new OptiCluster(data, params->metric, params->seed, cutoff, 0);
			method->SetClusterParameters(params->clusterParameters);
			ClusterExport* result = method->Execute();
			params->results.emplace_back(result);
			const ListVector& vector = result->GetListVector().listVector;
			SensSpecCalc calc(*data, &vector);
			calc.getResults(*data, tp, tn, fp, fn);
			params->fp += fp;
			params->fn += fn;
			params->tp += tp;
			params->tn += tn;
			delete method;
			delete data;
			continue;
		}
		RAbundVector rAbund = listVector.getRAbundVector();
		ClusterMethod* method = Utils::GetClusterMethod(params->clusterParameters.GetClusterType(), &listVector, &matrix,
			rAbund, params->seed, cutoff);
		method->SetClusterParameters(params->clusterParameters);
		params->results.emplace_back(method->Execute());
		delete method;

	}
}

// void ClusterSplit::runSensSpec() {
//
// }

void ClusterSplit::RunSenSpec() {
    std::vector<double> results;

    Sensitivity sens;
    double sensitivity = sens.getValue(tp, tn, fp, fn);
    results.push_back(sensitivity);
    Specificity spec;
    double specificity = spec.getValue(tp, tn, fp, fn);
    results.push_back(specificity);
    PPV ppv;
    double positivePredictiveValue = ppv.getValue(tp, tn, fp, fn);
    results.push_back(positivePredictiveValue);
    NPV npv;
    double negativePredictiveValue = npv.getValue(tp, tn, fp, fn);
    results.push_back(negativePredictiveValue);
    FDR fdr;
    double falseDiscoveryRate = fdr.getValue(tp, tn, fp, fn);
    results.push_back(falseDiscoveryRate);
    Accuracy acc;
    double accuracy = acc.getValue(tp, tn, fp, fn);
    results.push_back(accuracy);
    MCC mcc;
    double matthewsCorrCoef = mcc.getValue(tp, tn, fp, fn);
    results.push_back(matthewsCorrCoef);
    F1Score f1;
    double f1Score = f1.getValue(tp, tn, fp, fn);
    results.push_back(f1Score);

	std::vector<std::string> sensfileHeaders{"label","cutoff","tp","tn","fp","fn","sensitivity",
		"specificity","ppv","npv","fdr","accuracy","mcc","f1score"};

	std::string sensFile = std::to_string(cutoff) + ',' + std::to_string(cutoff) + ',' + std::to_string(tp) + ',' +
			std::to_string(tn) + ',' +
			std::to_string(fp) + ',' + std::to_string(fn) + ',';
	for (double result: results) { sensFile += std::to_string(result) + ','; }
	Utils::AddRowToDataFrameMap(dataframeMapSensMetrics, sensFile, sensfileHeaders);
}
