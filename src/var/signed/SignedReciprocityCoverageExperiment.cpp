/*
 * SignedReciprocityCoverageExperiment.cpp
 *
 *  Created on: Oct 14, 2010
 *      Author: duyvu
 */

#include "SignedReciprocityCoverageExperiment.h"

#include "EfficientSignedReciprocityParametricBoostrapSampler.h"

#include <iostream>
#include <fstream>
#include <string>

#define MATHLIB_STANDALONE
#include <Rmath.h>

namespace mixrerg {

SignedReciprocityCoverageExperiment::SignedReciprocityCoverageExperiment() {
	// TODO Auto-generated constructor stub

}

void SignedReciprocityCoverageExperiment::normalizeDyadProbs() {
	return;
}

SignedReciprocityCoverageExperiment::~SignedReciprocityCoverageExperiment() {
	// TODO Auto-generated destructor stub
	for (unsigned int k = 0; k < dyadProbs.size1(); k++)
		for (unsigned int l = 0; l < dyadProbs.size2(); l++)
			if (dyadProbs(k, l) != NULL)
				delete dyadProbs(k, l);
}

void SignedReciprocityCoverageExperiment::configure(
		const char* configurationFile, const char* delim) {

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

	if (configurationMap.count("Coverage_Sample_Size") == 1)
		coverageSampleSize = atoi(
				configurationMap["Coverage_Sample_Size"].c_str());
	else
		coverageSampleSize = 1000;

	if (configurationMap.count("EM_Configuration_File") == 1)
		emConfigurationFile = configurationMap["EM_Configuration_File"];
	else
		emConfigurationFile = "./";

	if (configurationMap.count("Bootstrap_Configuration_File") == 1)
		bootstrapConfigurationFile
				= configurationMap["Bootstrap_Configuration_File"];
	else
		bootstrapConfigurationFile = "./";

	if (configurationMap.count("Input_Directory") == 1)
		inputDirectory = configurationMap["Input_Directory"];
	else
		inputDirectory = "./";

	if (configurationMap.count("Output_Directory") == 1)
		outputDirectory = configurationMap["Output_Directory"];
	else
		outputDirectory = "./";

	if (configurationMap.count("Output_Precision") == 1)
		outputPrecision = atoi(configurationMap["Output_Precision"].c_str());
	else
		outputPrecision = 10;

	if (configurationMap.count("SEED_1") == 1)
		seed1 = atoi(configurationMap["SEED_1"].c_str());
	else
		seed1 = 5601874;

	if (configurationMap.count("SEED_2") == 1)
		seed2 = atoi(configurationMap["SEED_2"].c_str());
	else
		seed2 = 5710987;

	// configure EM Engine
	engine.configure(emConfigurationFile.c_str(), "\t");

	// configure Boostrap Sampler
	sampler = new EfficientSignedReciprocityParametricBoostrapSampler();
	sampler->configure(bootstrapConfigurationFile.c_str(), "\t");

	string inputAlphaFilename = inputDirectory + "/alpha.txt";
	ifstream inputAlphaFile;
	inputAlphaFile.open(inputAlphaFilename.c_str());
	if (inputAlphaFile.is_open()) {
		string line;
		// get the number of vertices and number of classes first
		getline(inputAlphaFile, line);
		StringTokenizer tokenizer = StringTokenizer(line, "\t");
		if (tokenizer.countTokens() == 2) {
			numOfVertices = tokenizer.nextIntToken();
			numOfClasses = tokenizer.nextIntToken();
			// get alpha values
			alpha = DoubleVector(numOfClasses);
			getline(inputAlphaFile, line);
			tokenizer = StringTokenizer(line, "\t");
			unsigned int index = 0;
			while (tokenizer.hasMoreTokens()) {
				alpha(index++) = tokenizer.nextDoubleToken();
			}
		} else
			cout << "There is no number of vertices or number of classes!!!"
					<< endl;
	} else
		cout << "Can not open the file " << inputAlphaFilename << endl;

	numOfPiParameterSets = 6;
	PIs = VectorOfDoubleMatrices(numOfPiParameterSets);
	for (unsigned int indexOfPI = 0; indexOfPI < numOfPiParameterSets; indexOfPI++) {
		PIs(indexOfPI) = DoubleMatrix(numOfClasses, numOfClasses);
		stringstream sIndexOfPI;
		sIndexOfPI << indexOfPI;
		string inputPIFilename = inputDirectory + "/PI" + sIndexOfPI.str()
				+ ".txt";
		sampler->readDoubleMatrixFromTextFile(PIs(indexOfPI),
				inputPIFilename.c_str());
	}

	// output CC
	DoubleVector overallClusteringCoefficientSample = DoubleVector(
			sampler->getNumberOfClusteringCoefficients());
	sampler->setParameters(numOfVertices, numOfClasses, alpha, PIs);
	sampler->computeOverallClusteringCoefficients(
			overallClusteringCoefficientSample, alpha, PIs);
	string outputOverallClusteringCoefficientSampleFilename = inputDirectory
			+ "/CC.txt";
	sampler->printDoubleVector2TextFile(overallClusteringCoefficientSample,
			outputOverallClusteringCoefficientSampleFilename.c_str(),
			outputPrecision);

	// Thetas
	VectorOfDoubleMatrices Thetas(sampler->getNumOfThetas());
	sampler->transformParameters(PIs, Thetas);
	for (unsigned int m = 0; m < Thetas.size(); m++) {
		stringstream indexOfTheta;
		indexOfTheta << m;
		string outputThetaFilename = inputDirectory + "/Theta"
				+ indexOfTheta.str() + ".txt";
		sampler->printDoubleMatrix2TextFile(Thetas(m),
				outputThetaFilename.c_str(), outputPrecision);
	}

	dyadProbs = DoublePointerMatrix(numOfClasses, numOfClasses);
	for (unsigned int k = 0; k < numOfClasses; k++)
		for (unsigned int l = 0; l < numOfClasses; l++) {
			dyadProbs(k, l) = new double[9];
			// y_ij = [-1,0,1,0,-1,1,-1,1,0];
			// y_ji = [0,-1,0,1,1,-1,-1,1,0];
			// p0(k,l),p0(l,k),p1(k,l),p1(l,k),p2(k,l),p2(l,k),p3(k,l),p4(k,l),p5(k,l)
			dyadProbs(k, l)[0] = PIs(0)(k, l);
			dyadProbs(k, l)[1] = PIs(0)(l, k);
			dyadProbs(k, l)[2] = PIs(1)(k, l);
			dyadProbs(k, l)[3] = PIs(1)(l, k);
			dyadProbs(k, l)[4] = PIs(2)(k, l);
			dyadProbs(k, l)[5] = PIs(2)(l, k);
			dyadProbs(k, l)[6] = PIs(3)(k, l);
			dyadProbs(k, l)[7] = PIs(4)(k, l);
			dyadProbs(k, l)[8] = PIs(5)(k, l);
			// should validate sum of these probabilities = 1 here
		}
	normalizeDyadProbs();
	cout << "Done with dyadProbs" << endl;
}

void SignedReciprocityCoverageExperiment::generateSample(
		NetworkSparseMatrix& Y, UnsignedIntVector& Z) {
	// copy alpha to an pointer-based array
	double * mixingProbs = new double[numOfClasses];
	for (unsigned int k = 0; k < alpha.size(); k++)
		mixingProbs[k] = alpha(k);

	// generate the latent class membership
	Z.clear();
	int * z_i = new int[numOfClasses];
	for (unsigned int i = 0; i < numOfVertices; i++) {
		rmultinom(1, mixingProbs, numOfClasses, z_i);
		int selectedClass = 0;
		while (z_i[selectedClass] == 0)
			selectedClass++;
		Z(i) = selectedClass;
	}

	// generate the adjacency matrix
	SparseGVOCoordinateV tempY = SparseGVOCoordinateV(numOfVertices,
			numOfVertices);
	int * dyad_ij = new int[9];
	for (unsigned int i = 0; i < numOfVertices; i++)
		for (unsigned int j = 0; j < i; j++) {
			rmultinom(1, dyadProbs(Z(i), Z(j)), 9, dyad_ij);
			int selectedDyadValue = 0;
			while (dyad_ij[selectedDyadValue] == 0)
				selectedDyadValue++;
			// y_ij = [-1,0,1,0,-1,1,-1,1,0];
			// y_ji = [0,-1,0,1,1,-1,-1,1,0];
			switch (selectedDyadValue) {
			// do not set 0 value because the matrix is in the sparse format
			// if we set value 0, it causes troubles when we iterate through nonzero elements
			// even being setting a value 0, an element turns to be nonzero ;-(
			case 0:
				tempY(i, j) = -1;
				break;
			case 1:
				tempY(j, i) = -1;
				break;
			case 2:
				tempY(i, j) = 1;
				break;
			case 3:
				tempY(j, i) = 1;
				break;
			case 4:
				tempY(i, j) = -1;
				tempY(j, i) = 1;
				break;
			case 5:
				tempY(i, j) = 1;
				tempY(j, i) = -1;
				break;
			case 6:
				tempY(i, j) = -1;
				tempY(j, i) = -1;
				break;
			case 7:
				tempY(i, j) = 1;
				tempY(j, i) = 1;
				break;
			case 8:
				break;
			default:
				cout << "ERROR!!! Undefined dyad value!" << endl;
			}
		}
	Y.assign(tempY);

	// clean up some pointer-based vectors and matrices
	if (mixingProbs != NULL)
		delete mixingProbs;
	if (z_i != NULL)
		delete z_i;
	if (dyad_ij != NULL)
		delete dyad_ij;
}

}
