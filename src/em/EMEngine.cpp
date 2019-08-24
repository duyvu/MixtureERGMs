/*
 * EMEngine.cpp
 *
 *  Created on: Jan 19, 2010
 *      Author: duyvu
 */

#include "EMEngine.h"

#include "model/binary/BinaryReciprocityModel.h"
#include "model/binary/MMBinaryReciprocityModel.h"

#include "model/signed/SignedReciprocityModel.h"
#include "model/signed/MMSignedReciprocityModel.h"

#include "util/StringTokenizer.h"

#include <map>

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include<ctime>
#include<time.h>
#include <math.h>

namespace mixrerg {

EMEngine::EMEngine() {
	// TODO Auto-generated constructor stub
}

EMEngine::~EMEngine() {
	// TODO Auto-generated destructor stub
	if (model != NULL)
		delete model;
}

void EMEngine::configure(const char* configurationFile, const char* delim) {

	std::cout << "The configuration file is " << configurationFile << std::endl;

	std::map<std::string, std::string> configurationMap;
	std::string line;
	std::ifstream fileStream(configurationFile);
	if (fileStream.is_open()) {
		while (!fileStream.eof()) {
			getline(fileStream, line);
			StringTokenizer strtok = StringTokenizer(line, delim);
			if (strtok.countTokens() >= 2) {
				std::string key = strtok.nextToken();
				std::string value = strtok.nextToken();
				configurationMap[key] = value;
			}
		}
	} else
		std::cout << "can not read the configuration file" << std::endl;

	for (std::map<std::string, std::string>::const_iterator it =
			configurationMap.begin(); it != configurationMap.end(); ++it) {
		std::cout << "(" << it->first << " , " << it->second << ")"
				<< std::endl;
	}

	int modelType = atoi(configurationMap["Model_Type"].c_str());
	switch (modelType) {
	case BINARY_GRAPH:
		model = new BinaryReciprocityModel();
		break;
	case MM_BINARY_GRAPH:
		model = new MMBinaryReciprocityModel();
		break;
	case SIGNED_GRAPH:
		model = new SignedReciprocityModel();
		break;
	case MM_SIGNED_GRAPH:
		model = new MMSignedReciprocityModel();
		break;
	default:
		std::cout << "Need to specify the model type. I am quitting now!"
				<< std::endl;
		return;
	}

	if (configurationMap.count("Data_File") == 1
			&& configurationMap.count("Data_Format") == 1) {
		const char* dataFile = configurationMap["Data_File"].c_str();
		int dataFormat = atoi(configurationMap["Data_Format"].c_str());
		std::string delim;
		switch (atoi(configurationMap["Delim"].c_str())) {
		case BLANK_SPACE:
			delim = " ";
			break;
		case TAB:
			delim = "\t";
			break;
		default:
			delim = " ";
			break;
		}
		model->setData(dataFile, dataFormat, delim.c_str());

		if (configurationMap.count("SEED_1") == 1)
			seed1 = atoi(configurationMap["SEED_1"].c_str());
		else
			seed1 = 45;

		if (configurationMap.count("SEED_2") == 1)
			seed2 = atoi(configurationMap["SEED_2"].c_str());
		else
			seed2 = 345;

		if (configurationMap.count("Num_Of_EM_Runs") == 1)
			numOfEMRuns = atoi(configurationMap["Num_Of_EM_Runs"].c_str());
		else
			numOfEMRuns = 1;

		if (configurationMap.count("Max_Num_Of_EM_Iterations") == 1)
			maxNumOfEMIterations = atoi(
					configurationMap["Max_Num_Of_EM_Iterations"].c_str());
		else
			maxNumOfEMIterations = 30;

		if (configurationMap.count("Max_Num_Of_E_Steps") == 1)
			maxNumOfESteps = atoi(
					configurationMap["Max_Num_Of_E_Steps"].c_str());
		else
			maxNumOfESteps = 30;

		if (configurationMap.count("Min_Tau") == 1)
			minTau = atof(configurationMap["Min_Tau"].c_str());
		else
			minTau = 1e-10;
		model->setMinTau(minTau);

		if (configurationMap.count("Min_Alpha") == 1)
			minAlpha = atof(configurationMap["Min_Alpha"].c_str());
		else
			minAlpha = 1e-6;
		model->setMinAlpha(minAlpha);

		if (configurationMap.count("Min_Pi") == 1)
			minPi = atof(configurationMap["Min_Pi"].c_str());
		else
			minPi = 1e-6;
		model->setMinPi(minPi);

		if (configurationMap.count("Tau_Precision") == 1)
			tauPrecision = atof(configurationMap["Tau_Precision"].c_str());
		else
			tauPrecision = 1e-6;

		if (configurationMap.count("ParamPrecision") == 1)
			paramPrecision = atof(configurationMap["ParamPrecision"].c_str());
		else
			paramPrecision = 1e-6;

		if (configurationMap.count("Min_Num_Of_Clusters") == 1)
			minNumOfClusters = atoi(
					configurationMap["Min_Num_Of_Clusters"].c_str());
		else
			minNumOfClusters = 3;

		if (configurationMap.count("Max_Num_Of_Clusters") == 1)
			maxNumOfClusters = atoi(
					configurationMap["Max_Num_Of_Clusters"].c_str());
		else
			maxNumOfClusters = 3;

	} else {
		std::cout << "No data available. I am quitting now!" << std::endl;
		return;
	}
}

void EMEngine::setData(const char* dataFile, int fileFormat,
		const char* delim) {
	model->setData(dataFile, fileFormat, delim);
}

void EMEngine::setData(const NetworkSparseMatrix& networkData) {
	model->setData(networkData);
}

void EMEngine::runEM(const int & _maxNumOfEMIterations,
		const int & _maxNumOfESteps, const double & _tauPrecision,
		const double & _paramPrecision) {
	std::cout << "EMEngine::runEM" << std::endl;

	cout.setf(ios::fixed, ios::floatfield);
	cout.precision(20);

	// assume that model has been created, the data and all parameters have been
	// configured
	model->resetParameters();
	double lastLLHLowerBound = 0;
	std::cout << "Starting a EM run" << std::endl;
	for (int emIteration = 0; emIteration < _maxNumOfEMIterations;
			emIteration++) {

		std::cout << "EM Iteration " << emIteration << std::endl;

		for (int eStep = 0; eStep < _maxNumOfESteps; eStep++) {
			model->runFixedPointEstimationEStep();
			if (!model->isTauSignificantlyChanged(_tauPrecision))
				break;
		}

		model->runModelEstimationMStep();
		if (!model->areAlphaPiSignificantlyChanged(_paramPrecision)) {
			cout << "alpha and pi are not significantly changed!" << endl;
			break;
		}

		// Covergence Details
		double newLLHLowerBound = model->getCompleteLogLikelihood();
		double largestTauChange = model->getLargestTauChange();

		// Convergence Criterion Check
		if ((emIteration > 0)
				&& (fabs(
						(newLLHLowerBound - lastLLHLowerBound)
								/ newLLHLowerBound) < REL_LLH))
			break;
		else
			lastLLHLowerBound = newLLHLowerBound;
	}
}

void EMEngine::runEMConvergenceDetails(const int & _maxNumOfEMIterations,
		const int & _maxNumOfESteps, const double & _tauPrecision,
		const double & _paramPrecision) {

	time_t startTime;
	time_t currentTime;

	std::cout << "EMEngine::runEM" << std::endl;

	cout.setf(ios::fixed, ios::floatfield);
	cout.precision(20);

	// assume that model has been created, the data and all parameters have been
	// configured
	model->resetParameters();

	double lastLLHLowerBound = model->getCompleteLogLikelihood();
	cout << "Iteration\t" << -1 << "\t" << lastLLHLowerBound << endl;

	currentTime = clock();
	startTime = currentTime;

	cout << "Time\t" << (double) (currentTime - startTime) / CLOCKS_PER_SEC << "\t"
			<< lastLLHLowerBound << endl;

	std::cout << "Starting a EM run" << std::endl;
	for (int emIteration = 0; emIteration < _maxNumOfEMIterations;
			emIteration++) {

		std::cout << "EM Iteration " << emIteration << std::endl;

		for (int eStep = 0; eStep < _maxNumOfESteps; eStep++) {
			std::cout << "E Step " << eStep << std::endl;
			model->runFixedPointEstimationEStep();
			if (!model->isTauSignificantlyChanged(_tauPrecision))
				break;
		}

		model->runModelEstimationMStep();
		if (!model->areAlphaPiSignificantlyChanged(_paramPrecision)) {
			cout << "alpha and pi are not significantly changed!" << endl;
			break;
		}

		// Covergence Details
		double newLLHLowerBound = model->getCompleteLogLikelihood();
		double largestTauChange = model->getLargestTauChange();
		cout << "Iteration\t" << emIteration << "\t" << newLLHLowerBound << "\t"
				<< largestTauChange << endl;
		currentTime = clock();
		cout << "Time\t" << (double) (currentTime - startTime) / CLOCKS_PER_SEC
				<< "\t" << newLLHLowerBound << endl;

		// Convergence Criterion Check
		if ((emIteration > 0)
				&& (fabs(
						(newLLHLowerBound - lastLLHLowerBound)
								/ newLLHLowerBound) < REL_LLH))
			break;
		else
			lastLLHLowerBound = newLLHLowerBound;
	}
}

void EMEngine::runEM() {
	std::cout << "EMEngine::runEM" << std::endl;

	cout.setf(ios::fixed, ios::floatfield);
	cout.precision(20);

	// assume that model has been created, the data and all parameters have been
	// configured
	model->resetParameters();

	std::cout << "Starting a EM run" << std::endl;
	for (int emIteration = 0; emIteration < maxNumOfEMIterations;
			emIteration++) {

		std::cout << "EM Iteration " << emIteration << std::endl;

		for (int eStep = 0; eStep < maxNumOfESteps; eStep++) {
			std::cout << "E Step " << eStep << std::endl;
			model->runFixedPointEstimationEStep();
			if (!model->isTauSignificantlyChanged(tauPrecision)) {
				cout << "tau is not significantly changed!" << endl;
				break;
			}
		}

		model->runModelEstimationMStep();
		if (!model->areAlphaPiSignificantlyChanged(paramPrecision)) {
			cout << "alpha and pi are not significantly changed!" << endl;
			break;
		}

		// Covergence Details
		double newLLHLowerBound = model->getCompleteLogLikelihood();
		double largestTauChange = model->getLargestTauChange();
		cout << emIteration << "\t";
		cout << newLLHLowerBound << "\t";
		cout << largestTauChange << endl;
	}
}

void EMEngine::runInitializedEM(int numOfClusters, const char* tauFile,
		DoubleMatrix& bestTau, DoubleVector& bestAlpha,
		VectorOfDoubleMatrices& bestPIs) {
	cout << "EMEngine::runInitializedEM" << endl;

	cout.setf(ios::fixed, ios::floatfield);
	cout.precision(20);

	// finalize all configuration setttings
	model->setModel(numOfClusters, minTau, minAlpha, minPi);

	// read initial tau
	UnsignedIntVector initTau = UnsignedIntVector(model->getNumOfVertices());
	ifstream tauStream(tauFile);
	if (tauStream.is_open()) {
		if (!tauStream.eof()) {
			string line;
			for (int i = 0; i < model->getNumOfVertices(); i++) {
				getline(tauStream, line);
				StringTokenizer strtok = StringTokenizer(line, "\t");
				initTau[i] = strtok.nextIntToken();
			}
			tauStream.close();
		}
	} else
		cout << "Unable to open tau file " << tauFile << endl;

	// initialize tau
	model->resetParameters(initTau);

	std::cout << "Starting a EM run" << std::endl;

	cout << "Init CLL = " << model->getCompleteLogLikelihood() << endl;

	for (int emIteration = 0; emIteration < maxNumOfEMIterations;
			emIteration++) {

		std::cout << "EM Iteration " << emIteration << std::endl;

		for (int eStep = 0; eStep < maxNumOfESteps; eStep++) {
			model->runFixedPointEstimationEStep();
			//			if (!model->isTauSignificantlyChanged(tauPrecision)) {
			//				cout << "tau is not significantly changed!" << endl;
			//				break;
			//			}
		}

		model->runModelEstimationMStep();
		//		if (!model->areAlphaPiSignificantlyChanged(paramPrecision)) {
		//			cout << "alpha and pi are not significantly changed!" << endl;
		//			break;
		//		}

		// Covergence Details
		double newLLHLowerBound = model->getCompleteLogLikelihood();
		double largestTauChange = model->getLargestTauChange();
		cout << emIteration << "\t";
		cout << newLLHLowerBound << "\t";
		cout << largestTauChange << endl;
	}

	model->copyParameters(bestTau, bestAlpha, bestPIs);
}

void EMEngine::getBestEMRun(const int _numOfClusters,
		const int & _maxNumOfEMIterations, const int & _maxNumOfESteps,
		const double & _tauPrecision, const double & _paramPrecision,
		const int & _numOfEMRuns, DoubleMatrix& bestTau,
		DoubleVector& bestAlpha, VectorOfDoubleMatrices& bestPIs,
		double& bestCLL, double& bestICL) {
	std::cout << "EMEngine::getBestEMRun" << std::endl;

	model->setModel(_numOfClusters, minTau, minAlpha, minPi);
	bestCLL = -DBL_MAX;
	for (int run = 0; run < _numOfEMRuns; run++) {
		std::cout << "Run " << run << ": " << std::endl;
		runEM(_maxNumOfEMIterations, _maxNumOfESteps, _tauPrecision,
				_paramPrecision);
		if (bestCLL < model->getCompleteLogLikelihood()) {
			std::cout << "Find the better" << std::endl;
			bestCLL = model->getCompleteLogLikelihood();
			std::cout << "better CLL = " << bestCLL << std::endl;
			bestICL = model->getICL();
			std::cout << "better ICL = " << bestICL << std::endl;
			model->copyParameters(bestTau, bestAlpha, bestPIs);
		}
	}
}

void EMEngine::getBestEMRunConvergenceDetails(const int _numOfClusters,
		const int & _maxNumOfEMIterations, const int & _maxNumOfESteps,
		const double & _tauPrecision, const double & _paramPrecision,
		const int & _numOfEMRuns, DoubleMatrix& bestTau,
		DoubleVector& bestAlpha, VectorOfDoubleMatrices& bestPIs,
		double& bestCLL, double& bestICL) {
	std::cout << "EMEngine::getBestEMRun" << std::endl;

	model->setModel(_numOfClusters, minTau, minAlpha, minPi);
	bestCLL = -DBL_MAX;
	for (int run = 0; run < _numOfEMRuns; run++) {
		std::cout << "Run " << run << ": " << std::endl;
		runEMConvergenceDetails(_maxNumOfEMIterations, _maxNumOfESteps,
				_tauPrecision, _paramPrecision);
		if (bestCLL < model->getCompleteLogLikelihood()) {
			std::cout << "Find the better" << std::endl;
			bestCLL = model->getCompleteLogLikelihood();
			std::cout << "better CLL = " << bestCLL << std::endl;
			bestICL = model->getICL();
			std::cout << "better ICL = " << bestICL << std::endl;
			model->copyParameters(bestTau, bestAlpha, bestPIs);
		}
	}
}

void EMEngine::getBestEMRun(const int _numOfClusters, DoubleMatrix& bestTau,
		DoubleVector& bestAlpha, VectorOfDoubleMatrices& bestPIs,
		double& bestCLL, double& bestICL) {
	std::cout << "EMEngine::getBestEMRun" << std::endl;

	model->setModel(_numOfClusters, minTau, minAlpha, minPi);
	bestCLL = -DBL_MAX;
	for (int run = 0; run < numOfEMRuns; run++) {
		std::cout << "Run " << run << ": " << std::endl;
		runEM();
		if (bestCLL < model->getCompleteLogLikelihood()) {
			std::cout << "Find the better" << std::endl;
			bestCLL = model->getCompleteLogLikelihood();
			std::cout << "better CLL = " << bestCLL << std::endl;
			bestICL = model->getICL();
			std::cout << "better ICL = " << bestICL << std::endl;
			model->copyParameters(bestTau, bestAlpha, bestPIs);
		}
	}
}

}
