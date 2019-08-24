/*
 * ParametricBoostrapSampler.cpp
 *
 *  Created on: Oct 4, 2010
 *      Author: duyvu
 */

#include "ParametricBoostrapSampler.h"

#include "model/binary/BinaryReciprocityModel.h"
#include "model/binary/MMBinaryReciprocityModel.h"
#include "model/binary/BootstrapBinaryReciprocityModel.h"

#include "model/signed/SignedReciprocityModel.h"
#include "model/signed/MMSignedReciprocityModel.h"
#include "model/signed/BootstrapSignedReciprocityModel.h"

#include "util/StringTokenizer.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cfloat>

namespace mixrerg {

ParametricBoostrapSampler::ParametricBoostrapSampler() {
	// TODO Auto-generated constructor stub
	model = NULL;
}

ParametricBoostrapSampler::~ParametricBoostrapSampler() {
	// TODO Auto-generated destructor stub
	if (model != NULL)
		delete model;
}

void ParametricBoostrapSampler::normalizeDyadProbs() {

}

void ParametricBoostrapSampler::generateSamples(unsigned int numOfSamples,
		const string& outputDirectory, unsigned outputPrecision) {

	// alpha
	ofstream outputAlphaFile;
	string outputAlphaFilename = outputDirectory + "/alpha.txt";
	outputAlphaFile.open(outputAlphaFilename.c_str());
	if (outputAlphaFile.is_open()) {
		outputAlphaFile.precision(outputPrecision);
		for (unsigned int k = 0; k < alpha.size(); k++)
			outputAlphaFile << alpha(k) << "\t";
		outputAlphaFile.close();
	} else
		cout << "Can not open the file " << outputAlphaFilename << endl;

	// PIs
	for (unsigned int m = 0; m < PIs.size(); m++) {
		ofstream outputPiFile;
		stringstream indexOfPi;
		indexOfPi << m;
		string outputPiFilename = outputDirectory + "/PI" + indexOfPi.str()
				+ ".txt";
		outputPiFile.open(outputPiFilename.c_str());
		if (outputPiFile.is_open()) {
			outputPiFile.precision(outputPrecision);
			for (unsigned int k = 0; k < numOfClasses; k++) {
				for (unsigned int l = 0; l < numOfClasses; l++)
					outputPiFile << PIs(m)(k, l) << "\t";
				outputPiFile << endl;
			}
			outputPiFile.close();
		} else
			cout << "Can not open the file " << outputPiFilename << endl;
	}

	// Thetas
	VectorOfDoubleMatrices Thetas(0);
	transformParameters(PIs, Thetas);
	for (unsigned int m = 0; m < Thetas.size(); m++) {
		stringstream indexOfTheta;
		indexOfTheta << m;
		string outputThetaFilename = outputDirectory + "/Theta"
				+ indexOfTheta.str() + ".txt";
		printDoubleMatrix2TextFile(Thetas(m), outputThetaFilename.c_str(),
				outputPrecision);
	}

	NetworkSparseMatrix Y = NetworkSparseMatrix(numOfVertices, numOfVertices);
	UnsignedIntVector Z = UnsignedIntVector(numOfVertices);
	for (unsigned int sampleIndex = 0; sampleIndex < numOfSamples; sampleIndex++) {

		std::cout << "Generating the sample " << sampleIndex << std::endl;

		// get the sample
		generateSample(Y, Z);

		// Z
		ofstream outputZFile;
		stringstream sSampleIndex;
		sSampleIndex << sampleIndex;
		string outputZFilename = outputDirectory + "/Z_" + sSampleIndex.str()
				+ ".txt";
		outputZFile.open(outputZFilename.c_str());
		if (outputZFile.is_open()) {
			outputZFile << numOfVertices << "\t" << numOfClasses << endl;
			for (unsigned int i = 0; i < numOfVertices; i++)
				outputZFile << Z(i) << endl;
			outputZFile.close();
		} else
			cout << "Can not open the file " << outputZFilename << endl;

		// Y
		ofstream outputYFile;
		string outputYFilename = outputDirectory + "/Y_" + sSampleIndex.str()
				+ ".txt";
		outputYFile.open(outputYFilename.c_str());
		if (outputYFile.is_open()) {
			outputYFile << numOfVertices << endl;
			for (NetworkSparseMatrix::iterator1 i1 = Y.begin1(); i1 != Y.end1(); ++i1) {
				for (NetworkSparseMatrix::iterator2 i2 = i1.begin(); i2
						!= i1.end(); ++i2) {
					unsigned int i = i2.index1();
					unsigned int j = i2.index2();
					outputYFile << i << "\t" << j << "\t" << Y(i, j) << endl;
				}
			}
			outputYFile.close();
		} else
			cout << "Can not open the file " << outputYFilename << endl;
	}
}

void ParametricBoostrapSampler::generateSamples(unsigned fromSampleIndex,
		unsigned toSampleIndex, const string& outputDirectory,
		unsigned outputPrecision) {

	// alpha
	ofstream outputAlphaFile;
	string outputAlphaFilename = outputDirectory + "/alpha.txt";
	outputAlphaFile.open(outputAlphaFilename.c_str());
	if (outputAlphaFile.is_open()) {
		outputAlphaFile.precision(outputPrecision);
		for (unsigned int k = 0; k < alpha.size(); k++)
			outputAlphaFile << alpha(k) << "\t";
		outputAlphaFile.close();
	} else
		cout << "Can not open the file " << outputAlphaFilename << endl;

	// PIs
	for (unsigned int m = 0; m < PIs.size(); m++) {
		ofstream outputPiFile;
		stringstream indexOfPi;
		indexOfPi << m;
		string outputPiFilename = outputDirectory + "/PI" + indexOfPi.str()
				+ ".txt";
		outputPiFile.open(outputPiFilename.c_str());
		if (outputPiFile.is_open()) {
			outputPiFile.precision(outputPrecision);
			for (unsigned int k = 0; k < numOfClasses; k++) {
				for (unsigned int l = 0; l < numOfClasses; l++)
					outputPiFile << PIs(m)(k, l) << "\t";
				outputPiFile << endl;
			}
			outputPiFile.close();
		} else
			cout << "Can not open the file " << outputPiFilename << endl;
	}

	// Thetas
	VectorOfDoubleMatrices Thetas(0);
	transformParameters(PIs, Thetas);
	for (unsigned int m = 0; m < Thetas.size(); m++) {
		stringstream indexOfTheta;
		indexOfTheta << m;
		string outputThetaFilename = outputDirectory + "/Theta"
				+ indexOfTheta.str() + ".txt";
		printDoubleMatrix2TextFile(Thetas(m), outputThetaFilename.c_str(),
				outputPrecision);
	}

	NetworkSparseMatrix Y = NetworkSparseMatrix(numOfVertices, numOfVertices);
	UnsignedIntVector Z = UnsignedIntVector(numOfVertices);
	for (unsigned int sampleIndex = fromSampleIndex; sampleIndex
			<= toSampleIndex; sampleIndex++) {

		std::cout << "Generating the sample " << sampleIndex << std::endl;

		// get the sample
		generateSample(Y, Z);

		// Z
		ofstream outputZFile;
		stringstream sSampleIndex;
		sSampleIndex << sampleIndex;
		string outputZFilename = outputDirectory + "/Z_" + sSampleIndex.str()
				+ ".txt";
		outputZFile.open(outputZFilename.c_str());
		if (outputZFile.is_open()) {
			outputZFile << numOfVertices << "\t" << numOfClasses << endl;
			for (unsigned int i = 0; i < numOfVertices; i++)
				outputZFile << Z(i) << endl;
			outputZFile.close();
		} else
			cout << "Can not open the file " << outputZFilename << endl;

		// Y
		ofstream outputYFile;
		string outputYFilename = outputDirectory + "/Y_" + sSampleIndex.str()
				+ ".txt";
		outputYFile.open(outputYFilename.c_str());
		if (outputYFile.is_open()) {
			outputYFile << numOfVertices << endl;
			for (NetworkSparseMatrix::iterator1 i1 = Y.begin1(); i1 != Y.end1(); ++i1) {
				for (NetworkSparseMatrix::iterator2 i2 = i1.begin(); i2
						!= i1.end(); ++i2) {
					unsigned int i = i2.index1();
					unsigned int j = i2.index2();
					outputYFile << i << "\t" << j << "\t" << Y(i, j) << endl;
				}
			}
			outputYFile.close();
		} else
			cout << "Can not open the file " << outputYFilename << endl;
	}
}

void ParametricBoostrapSampler::configure(const char* configurationFile,
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

	modelType = atoi(configurationMap["Model_Type"].c_str());
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
	case BOOTSTRAP_BINARY_GRAPH:
		std::cout << "Using BOOTSTRAP_BINARY_GRAPH!" << std::endl;
		model = new BootstrapBinaryReciprocityModel();
		break;
	case BOOTSTRAP_SIGNED_GRAPH:
		std::cout << "Using BOOTSTRAP_SIGNED_GRAPH!" << std::endl;
		model = new BootstrapSignedReciprocityModel();
		break;
	default:
		std::cout << "Model type must be specified!!!" << std::endl;
		exit(-1);
	}

	if (configurationMap.count("SEED_1") == 1)
		seed1 = atoi(configurationMap["SEED_1"].c_str());
	else
		seed1 = 5601874;

	if (configurationMap.count("SEED_2") == 1)
		seed2 = atoi(configurationMap["SEED_2"].c_str());
	else
		seed2 = 5710987;

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

	if (configurationMap.count("Bootstrap_Sample_Size") == 1)
		bootstrapSampleSize = atoi(
				configurationMap["Bootstrap_Sample_Size"].c_str());
	else
		bootstrapSampleSize = 1000;

	if (configurationMap.count("Num_Of_EM_Runs") == 1)
		numOfEMRuns = atoi(configurationMap["Num_Of_EM_Runs"].c_str());
	else
		numOfEMRuns = 1;

	if (configurationMap.count("Max_Num_Of_EM_Iterations") == 1)
		maxNumOfEMIterations = atoi(
				configurationMap["Max_Num_Of_EM_Iterations"].c_str());
	else
		maxNumOfEMIterations = 30;

	// if modelType is MM_SIGNED_GRAPH
	if (modelType == MM_SIGNED_GRAPH)
		maxNumOfESteps = 1;
	else {
		if (configurationMap.count("Max_Num_Of_E_Steps") == 1)
			maxNumOfESteps = atoi(
					configurationMap["Max_Num_Of_E_Steps"].c_str());
		else
			maxNumOfESteps = 30;
	}

	if (configurationMap.count("Tau_Precision") == 1)
		tauPrecision = atof(configurationMap["Tau_Precision"].c_str());
	else
		tauPrecision = 1e-6;

	if (configurationMap.count("ParamPrecision") == 1)
		paramPrecision = atof(configurationMap["ParamPrecision"].c_str());
	else
		paramPrecision = 1e-6;
}

void ParametricBoostrapSampler::runEM(const UnsignedIntVector& Z) {
	std::cout << "ParametricBoostrapSampler::runEM" << std::endl;

	// assume that model has been created
	// i.e. the data and all parameters have been configured
	model->resetParameters(Z);
	double initCLL = model->getCompleteLogLikelihood();
	cout << "initCLL = " << initCLL << endl;

	std::cout << "Starting a EM run" << std::endl;
	for (int emIteration = 0; emIteration < maxNumOfEMIterations; emIteration++) {

		std::cout << "EM Iteration " << emIteration << std::endl;

		for (int eStep = 0; eStep < maxNumOfESteps; eStep++) {
			std::cout << "E Step " << eStep << std::endl;
			model->runFixedPointEstimationEStep();
			if (!model->isTauSignificantlyChanged(tauPrecision))
				break;
		}

		model->runModelEstimationMStep();
		if (!model->areAlphaPiSignificantlyChanged(paramPrecision)) {
			break;
		}

		double currentCLL = model->getCompleteLogLikelihood();
		cout << "currentCLL = " << currentCLL << endl;

		if ((fabs(initCLL - currentCLL) / currentCLL) < 1e-10)
			break;
	}
}

void ParametricBoostrapSampler::getBestEMRun(const UnsignedIntVector& Z,
		DoubleMatrix& bestTau, DoubleVector& bestAlpha,
		VectorOfDoubleMatrices& bestPIs, double& bestCLL, double& bestICL) {
	std::cout << "ParametricBoostrapSampler::getBestEMRun" << std::endl;
	bestCLL = -DBL_MAX;
	for (int run = 0; run < numOfEMRuns; run++) {
		runEM(Z);
		if (numOfEMRuns <= 1) {
			model->copyParameters(bestTau, bestAlpha, bestPIs);
		} else if (bestCLL < model->getCompleteLogLikelihood()) {
			std::cout << "Find a better model!!!" << std::endl;
			bestCLL = model->getCompleteLogLikelihood();
			std::cout << "better CLL = " << bestCLL << std::endl;
			bestICL = model->getICL();
			std::cout << "better ICL = " << bestICL << std::endl;
			model->copyParameters(bestTau, bestAlpha, bestPIs);
		}
	}
}

void ParametricBoostrapSampler::estimateSample(unsigned sampleIndex,
		const string& inputDirectory, const string& outputDirectory,
		unsigned int outputPrecision) {

	clock_t start = clock();

	stringstream sSampleIndex;
	sSampleIndex << sampleIndex;

	string inputYFilename = inputDirectory + "/Y_" + sSampleIndex.str()
			+ ".txt";
	model->setData(inputYFilename.c_str(), 0, "\t");

	model->setModel(numOfClasses, minTau, minAlpha, minPi);

	UnsignedIntVector Z = UnsignedIntVector(numOfVertices);
	ifstream inputZFile;
	string inputZFilename = inputDirectory + "/Z_" + sSampleIndex.str()
			+ ".txt";
	cout << "Z file is " << inputZFilename << endl;
	inputZFile.open(inputZFilename.c_str());
	if (inputZFile.is_open()) {
		string line;
		getline(inputZFile, line); // just ignore the numbers of vertices and classes
		unsigned int index = 0;
		while (!inputZFile.eof() && index < numOfVertices) {
			getline(inputZFile, line);
			stringstream ss(line);
			ss >> Z(index);
			index++;
		}
		if (index != numOfVertices)
			cout << "There are some missing membership values!!!" << endl;
		inputZFile.close();
	} else
		cout << "Can not open the file " << inputZFilename << endl;

	DoubleMatrix bestTau = DoubleMatrix(numOfVertices, numOfClasses);
	DoubleVector bestAlpha = DoubleVector(numOfClasses);
	VectorOfDoubleMatrices bestPIs = VectorOfDoubleMatrices(
			numOfPiParameterSets);
	for (unsigned int m = 0; m < numOfPiParameterSets; m++)
		bestPIs(m) = DoubleMatrix(numOfClasses, numOfClasses);
	double bestCLL = 0;
	double bestICL = 0;
	getBestEMRun(Z, bestTau, bestAlpha, bestPIs, bestCLL, bestICL);

	// IF WE NEED TO RESOLVE LABEL SWITCHING, DO IT HERE

	// print out bestTau, bestAlpha, bestPIs, bestCLL, bestICL into files

	// alpha
	ofstream outputAlphaFile;
	string outputAlphaFilename = outputDirectory + "/alpha_"
			+ sSampleIndex.str() + ".txt";
	outputAlphaFile.open(outputAlphaFilename.c_str());
	if (outputAlphaFile.is_open()) {
		outputAlphaFile.precision(outputPrecision);
		outputAlphaFile << bestCLL << "\t" << bestICL << endl;
		printDoubleVector2TextFile(bestAlpha, outputAlphaFile);
		outputAlphaFile.close();
	} else
		cout << "Can not open the file " << outputAlphaFilename << endl;

	// PIs
	for (unsigned int m = 0; m < bestPIs.size(); m++) {
		stringstream indexOfPi;
		indexOfPi << m;
		string outputPiFilename = outputDirectory + "/PI" + indexOfPi.str()
				+ "_" + sSampleIndex.str() + ".txt";
		printDoubleMatrix2TextFile(bestPIs(m), outputPiFilename.c_str(),
				outputPrecision);
	}

	// Thetas
	VectorOfDoubleMatrices bestThetas(0);
	transformParameters(bestPIs, bestThetas);
	for (unsigned int m = 0; m < bestThetas.size(); m++) {
		stringstream indexOfTheta;
		indexOfTheta << m;
		string outputThetaFilename = outputDirectory + "/Theta"
				+ indexOfTheta.str() + "_" + sSampleIndex.str() + ".txt";
		printDoubleMatrix2TextFile(bestThetas(m), outputThetaFilename.c_str(),
				outputPrecision);
	}

	// tau
	string outputTauFilename = outputDirectory + "/tau_" + sSampleIndex.str()
			+ ".txt";
	printDoubleMatrix2TextFile(bestTau, outputTauFilename.c_str(),
			outputPrecision);

	clock_t jobTime = clock() - start;
	ofstream outputJobTimeFile;
	string outputJobTimeFilename = outputDirectory + "/JobTimes.txt";
	outputJobTimeFile.open(outputJobTimeFilename.c_str(), ios::app);
	if (outputJobTimeFile.is_open()) {
		time_t rawtime;
		time(&rawtime);
		jobTime = jobTime / (1e6); // to seconds
		unsigned int hours = jobTime / (60 * 60);
		jobTime = jobTime % (60 * 60);
		unsigned int minutes = jobTime / 60;
		unsigned int seconds = jobTime % 60;
		stringstream sJobTime;
		sJobTime << hours << ":";
		sJobTime << minutes << ":";
		sJobTime << seconds;
		// Sample index + Job time in hours:minutes:seconds s+ Done time
		outputJobTimeFile << sampleIndex << "\t";
		outputJobTimeFile << sJobTime.str() << "\t";
		outputJobTimeFile << ctime(&rawtime);
		outputJobTimeFile.close();
	} else
		cout << "Can not open the file " << outputJobTimeFilename << endl;
}

void ParametricBoostrapSampler::estimateSamples(unsigned fromSampleIndex,
		unsigned toSampleIndex, const string& inputDirectory,
		const string& outputDirectory, unsigned int outputPrecision) {
	for (unsigned sampleIndex = fromSampleIndex; sampleIndex <= toSampleIndex; sampleIndex++)
		estimateSample(sampleIndex, inputDirectory, outputDirectory,
				outputPrecision);
}

void ParametricBoostrapSampler::summarizeThetaSamples(unsigned fromSampleIndex,
		unsigned toSampleIndex, const string& inputDirectory,
		const string& outputDirectory, unsigned int outputPrecision) {

	unsigned int numOfThetas = getNumOfThetas();
	VectorOfMatrixOfDoubleVectors results = VectorOfMatrixOfDoubleVectors(
			numOfThetas);
	unsigned int sampleSize = toSampleIndex - fromSampleIndex + 1;
	for (unsigned int m = 0; m < numOfThetas; m++) {
		results(m) = MatrixOfDoubleVectors(numOfClasses, numOfClasses);
		for (unsigned int k = 0; k < numOfClasses; k++)
			for (unsigned int l = 0; l < numOfClasses; l++)
				results(m)(k, l) = DoubleVector(sampleSize);
	}

	DoubleMatrix matrix = DoubleMatrix(numOfClasses, numOfClasses);
	for (unsigned int sampleIndex = fromSampleIndex, i = 0; sampleIndex
			<= toSampleIndex; sampleIndex++, i++) {
		stringstream sSampleIndex;
		sSampleIndex << sampleIndex;
		for (unsigned int m = 0; m < numOfThetas; m++) {
			stringstream indexOfTheta;
			indexOfTheta << m;
			string inputThetaFilename = inputDirectory + "/Theta"
					+ indexOfTheta.str() + "_" + sSampleIndex.str() + ".txt";
			readDoubleMatrixFromTextFile(matrix, inputThetaFilename.c_str());
			for (unsigned int k = 0; k < numOfClasses; k++)
				for (unsigned int l = 0; l < numOfClasses; l++)
					results(m)(k, l)(i) = matrix(k, l);
		}
	}

	for (unsigned int m = 0; m < numOfThetas; m++) {
		stringstream indexOfTheta;
		indexOfTheta << m;
		for (unsigned int k = 0; k < numOfClasses; k++) {
			stringstream sK;
			sK << k;
			for (unsigned int l = 0; l < numOfClasses; l++) {
				stringstream sL;
				sL << l;
				string outputThetaFilename = outputDirectory + "/Theta"
						+ indexOfTheta.str() + "_c" + sK.str() + "_c"
						+ sL.str() + ".txt";
				cout << "Output the file " << outputThetaFilename << endl;
				printDoubleVector2TextFile(results(m)(k, l),
						outputThetaFilename.c_str(), outputPrecision);
			}
		}
	}
}

void ParametricBoostrapSampler::obtainThetaBootstrapSamples(
		const string& outputDirectory, const unsigned int outputPrecision) {

	unsigned int numOfThetas = getNumOfThetas();
	VectorOfMatrixOfDoubleVectors bootstrapThetaSamples =
			VectorOfMatrixOfDoubleVectors(numOfThetas);
	for (unsigned int m = 0; m < numOfThetas; m++) {
		bootstrapThetaSamples(m) = MatrixOfDoubleVectors(numOfClasses,
				numOfClasses);
		for (unsigned int k = 0; k < numOfClasses; k++)
			for (unsigned int l = 0; l < numOfClasses; l++)
				bootstrapThetaSamples(m)(k, l) = DoubleVector(
						bootstrapSampleSize);
	}

	DoubleMatrix bootstrapOverallClusteringCoefficientSamples = DoubleMatrix(
			bootstrapSampleSize, getNumberOfClusteringCoefficients());

	// output alpha, PIs, and Thetas

	// alpha
	ofstream outputAlphaFile;
	string outputAlphaFilename = outputDirectory + "/alpha.txt";
	outputAlphaFile.open(outputAlphaFilename.c_str());
	if (outputAlphaFile.is_open()) {
		outputAlphaFile.precision(outputPrecision);
		for (unsigned int k = 0; k < alpha.size(); k++)
			outputAlphaFile << alpha(k) << "\t";
		outputAlphaFile.close();
	} else
		cout << "Can not open the file " << outputAlphaFilename << endl;

	// PIs
	for (unsigned int m = 0; m < PIs.size(); m++) {
		ofstream outputPiFile;
		stringstream indexOfPi;
		indexOfPi << m;
		string outputPiFilename = outputDirectory + "/PI" + indexOfPi.str()
				+ ".txt";
		outputPiFile.open(outputPiFilename.c_str());
		if (outputPiFile.is_open()) {
			outputPiFile.precision(outputPrecision);
			for (unsigned int k = 0; k < numOfClasses; k++) {
				for (unsigned int l = 0; l < numOfClasses; l++)
					outputPiFile << PIs(m)(k, l) << "\t";
				outputPiFile << endl;
			}
			outputPiFile.close();
		} else
			cout << "Can not open the file " << outputPiFilename << endl;
	}

	// output CC
	DoubleVector overallClusteringCoefficientSample = DoubleVector(
			getNumberOfClusteringCoefficients());
	computeOverallClusteringCoefficients(overallClusteringCoefficientSample,
			alpha, PIs);
	string outputOverallClusteringCoefficientSampleFilename = outputDirectory
			+ "/CC_Original.txt";
	printDoubleVector2TextFile(overallClusteringCoefficientSample,
			outputOverallClusteringCoefficientSampleFilename.c_str(),
			outputPrecision);

	// Thetas
	VectorOfDoubleMatrices Thetas(numOfThetas);
	transformParameters(PIs, Thetas);
	for (unsigned int m = 0; m < Thetas.size(); m++) {
		stringstream indexOfTheta;
		indexOfTheta << m;
		string outputThetaFilename = outputDirectory + "/Theta"
				+ indexOfTheta.str() + ".txt";
		printDoubleMatrix2TextFile(Thetas(m), outputThetaFilename.c_str(),
				outputPrecision);
	}

	NetworkSparseMatrix Y = NetworkSparseMatrix(numOfVertices, numOfVertices);
	UnsignedIntVector Z = UnsignedIntVector(numOfVertices);
	DoubleMatrix bootstrapTau = DoubleMatrix(numOfVertices, numOfClasses);
	DoubleVector bootstrapAlpha = DoubleVector(numOfClasses);
	VectorOfDoubleMatrices bootstrapPIs = VectorOfDoubleMatrices(
			numOfPiParameterSets);
	for (unsigned int m = 0; m < numOfPiParameterSets; m++)
		bootstrapPIs(m) = DoubleMatrix(numOfClasses, numOfClasses);
	DoubleVector bootstrapOverallClusteringCoefficientSample = DoubleVector(
			getNumberOfClusteringCoefficients());
	double bootstrapCLL = 0;
	double bootstrapICL = 0;
	for (unsigned int sampleIndex = 0; sampleIndex < bootstrapSampleSize; sampleIndex++) {

		// get the sample
		std::cout << "Generating the sample " << sampleIndex << std::endl;
		generateSample(Y, Z);

		// reset data
		cout << "Resetting data" << endl;
		model->setData(Y);

		// reset model
		cout << "Resetting model" << endl;
		model->setModel(numOfClasses, minTau, minAlpha, minPi);

		// estimate parameters from the sample

		getBestEMRun(Z, bootstrapTau, bootstrapAlpha, bootstrapPIs,
				bootstrapCLL, bootstrapICL);

		//printAlpha(bootstrapAlpha);
		//printPIs(bootstrapPIs);

		// transform PIs into Thetas
		VectorOfDoubleMatrices bootstrapThetas(numOfThetas);
		transformParameters(bootstrapPIs, bootstrapThetas);
		for (unsigned int m = 0; m < bootstrapThetas.size(); m++)
			for (unsigned int k = 0; k < numOfClasses; k++)
				for (unsigned int l = 0; l < numOfClasses; l++)
					bootstrapThetaSamples(m)(k, l)(sampleIndex)
							= bootstrapThetas(m)(k, l);

		computeOverallClusteringCoefficients(
				bootstrapOverallClusteringCoefficientSample, bootstrapAlpha,
				bootstrapPIs);
		for (unsigned int k = 0; k < getNumberOfClusteringCoefficients(); k++)
			bootstrapOverallClusteringCoefficientSamples(sampleIndex, k)
					= bootstrapOverallClusteringCoefficientSample(k);

		std::cout << "Done with the sample " << sampleIndex << std::endl;
	}

	for (unsigned int m = 0; m < numOfThetas; m++) {
		stringstream indexOfTheta;
		indexOfTheta << m;
		for (unsigned int k = 0; k < numOfClasses; k++) {
			stringstream sK;
			sK << k;
			for (unsigned int l = 0; l < numOfClasses; l++) {
				stringstream sL;
				sL << l;
				string outputThetaFilename = outputDirectory + "/Theta"
						+ indexOfTheta.str() + "_c" + sK.str() + "_c"
						+ sL.str() + ".txt";
				cout << "Output the file " << outputThetaFilename << endl;
				printDoubleVector2TextFile(bootstrapThetaSamples(m)(k, l),
						outputThetaFilename.c_str(), outputPrecision);
			}
		}
	}

	string outputClusteringCoefficientFilename = outputDirectory + "/CC.txt";
	cout << "Output the file " << outputClusteringCoefficientFilename << endl;
	printDoubleMatrix2TextFile(bootstrapOverallClusteringCoefficientSamples,
			outputClusteringCoefficientFilename.c_str(), outputPrecision);
}

void ParametricBoostrapSampler::printAlpha(const DoubleVector& alpha) {
	cout.precision(6);
	cout << "Printing ALPHA: " << endl;
	cout << alpha(0);
	unsigned int numOfClasses = alpha.size();
	for (unsigned int k = 1; k < numOfClasses; k++)
		cout << ",  " << alpha(k);
	cout << endl;
}

void ParametricBoostrapSampler::printPIs(const VectorOfDoubleMatrices& PIs) {
	cout.precision(6);
	unsigned int numOfClasses = PIs(0).size1();
	if (modelType == BINARY_GRAPH || modelType == MM_BINARY_GRAPH || modelType
			== BOOTSTRAP_BINARY_GRAPH) {
		cout << "pi10: " << endl;
		for (unsigned int k = 0; k < numOfClasses; k++) {
			for (unsigned int l = 0; l < numOfClasses; l++)
				cout << PIs(0)(k, l) << "  ";
			cout << endl;
		}
		cout << "pi11: " << endl;
		for (unsigned int k = 0; k < numOfClasses; k++) {
			for (unsigned int l = 0; l < numOfClasses; l++)
				cout << PIs(1)(k, l) << "  ";
			cout << endl;
		}
		cout << "pi00: " << endl;
		for (unsigned int k = 0; k < numOfClasses; k++) {
			for (unsigned int l = 0; l < numOfClasses; l++)
				cout << PIs(2)(k, l) << "  ";
			cout << endl;
		}
	} else if (modelType == SIGNED_GRAPH || modelType == MM_SIGNED_GRAPH
			|| modelType == BOOTSTRAP_SIGNED_GRAPH) {
		cout << "pi1: " << endl;
		for (unsigned int k = 0; k < numOfClasses; k++) {
			for (unsigned int l = 0; l < numOfClasses; l++)
				cout << PIs(0)(k, l) << "  ";
			cout << endl;
		}
		cout << "pi2: " << endl;
		for (unsigned int k = 0; k < numOfClasses; k++) {
			for (unsigned int l = 0; l < numOfClasses; l++)
				cout << PIs(1)(k, l) << "  ";
			cout << endl;
		}
		cout << "pi3: " << endl;
		for (unsigned int k = 0; k < numOfClasses; k++) {
			for (unsigned int l = 0; l < numOfClasses; l++)
				cout << PIs(2)(k, l) << "  ";
			cout << endl;
		}
		cout << "pi4m: " << endl;
		for (unsigned int k = 0; k < numOfClasses; k++) {
			for (unsigned int l = 0; l < numOfClasses; l++)
				cout << PIs(3)(k, l) << "  ";
			cout << endl;
		}
		cout << "pi4p: " << endl;
		for (unsigned int k = 0; k < numOfClasses; k++) {
			for (unsigned int l = 0; l < numOfClasses; l++)
				cout << PIs(4)(k, l) << "  ";
			cout << endl;
		}
		cout << "pi5: " << endl;
		for (unsigned int k = 0; k < numOfClasses; k++) {
			for (unsigned int l = 0; l < numOfClasses; l++)
				cout << PIs(5)(k, l) << "  ";
			cout << endl;
		}
	}
}

}
