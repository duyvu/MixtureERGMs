/*
 * DelayedDataEMEngine.cpp
 *
 *  Created on: Oct 14, 2010
 *      Author: duyvu
 */

#include "DelayedDataEMEngine.h"

#include "model/binary/BinaryReciprocityModel.h"
#include "model/binary/MMBinaryReciprocityModel.h"

#include "model/signed/SignedReciprocityModel.h"
#include "model/signed/MMSignedReciprocityModel.h"

#include "util/StringTokenizer.h"

#include <iostream>
#include <fstream>

namespace mixrerg {

DelayedDataEMEngine::DelayedDataEMEngine() {
	// TODO Auto-generated constructor stub

}

DelayedDataEMEngine::~DelayedDataEMEngine() {
	// TODO Auto-generated destructor stub
}

void DelayedDataEMEngine::configure(const char* configurationFile,
		const char* delim) {

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
		maxNumOfESteps = atoi(configurationMap["Max_Num_Of_E_Steps"].c_str());
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
		minNumOfClusters
				= atoi(configurationMap["Min_Num_Of_Clusters"].c_str());
	else
		minNumOfClusters = 3;

	if (configurationMap.count("Max_Num_Of_Clusters") == 1)
		maxNumOfClusters
				= atoi(configurationMap["Max_Num_Of_Clusters"].c_str());
	else
		maxNumOfClusters = 3;
}

}
