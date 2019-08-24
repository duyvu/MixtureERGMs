//============================================================================
// Name        : mixrerg.cpp
// Author      : Duy Vu
// Version     :
// Copyright   : Penn State Copyright
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "utils.h"

#include "model/DataTypes.h"

#include "model/ReciprocityModel.h"

#include "model/binary/BinaryReciprocityModel.h"
#include "model/binary/MMBinaryReciprocityModel.h"

#include "model/signed/SignedReciprocityModel.h"
#include "model/signed/MMSignedReciprocityModel.h"

#include "em/EMEngine.h"

#include "var/ParametricBoostrapSampler.h"

#include "var/binary/BinaryReciprocityParametricBoostrapSampler.h"
#include "var/binary/EfficientBinaryReciprocityParametricBoostrapSampler.h"

#include "var/signed/SignedReciprocityParametricBoostrapSampler.h"
#include "var/signed/EfficientSignedReciprocityParametricBoostrapSampler.h"

#include "var/signed/SignedReciprocityCoverageExperiment.h"
#include "var/signed/EfficientSignedReciprocityCoverageExperiment.h"

#define MATHLIB_STANDALONE
#include <Rmath.h>

#include <iostream>
#include <fstream>
#include <ctime>
#include <string>

using namespace std;
using namespace mixrerg;

void estimateParameters(const char* confFile, const char* outputDir,
		unsigned int numOfClusters, unsigned int fromRun, unsigned int toRun) {

	cout << "Testing Signed Reciprocity Model " << endl;
	cout << "confFile is " << confFile << endl;
	cout << "outputDirectory is " << outputDir << endl;
	cout << "numOfClusters is " << numOfClusters << endl;

	EMEngine engine = EMEngine();
	engine.configure(confFile, "\t");
	engine.setMinNumOfClusters(numOfClusters);
	engine.setMaxNumOfClusters(numOfClusters);
	engine.setNumOfEMRuns(toRun - fromRun + 1);
	set_seed(fromRun, toRun);
	cout << "Finished configuration" << endl;

	DoubleMatrix bestTau = DoubleMatrix(engine.getNumOfVertices(),
			engine.getMinNumOfClusters());
	DoubleVector bestAlpha = DoubleVector(engine.getMinNumOfClusters());
	VectorOfDoubleMatrices bestPIs = VectorOfDoubleMatrices(
			engine.getNumOfParameterSets());
	for (int t = 0; t < engine.getNumOfParameterSets(); t++) {
		bestPIs(t) = DoubleMatrix(engine.getMinNumOfClusters(),
				engine.getMinNumOfClusters());
	}
	double bestCLL = 0;
	double bestICL = 0;

	engine.getBestEMRun(engine.getMinNumOfClusters(),
			engine.getMaxNumOfEMIterations(), engine.getMaxNumOfESteps(),
			engine.getTauPrecision(), engine.getParamPrecision(),
			engine.getNumOfEMRuns(), bestTau, bestAlpha, bestPIs, bestCLL,
			bestICL);

	string outputDirectory = string(outputDir);
	unsigned int outputPrecision = 10;

	stringstream postfix;
	postfix << engine.getMinNumOfClusters() << "_";
	postfix << fromRun << "_";
	postfix << toRun << "_";
	postfix << engine.getMaxNumOfEMIterations() << "_";
	postfix << engine.getMaxNumOfESteps();

	string outputAlphaFilename = outputDirectory + "/alpha_" + postfix.str()
			+ ".txt";
	ofstream outputAlphaFile;
	outputAlphaFile.open(outputAlphaFilename.c_str());
	if (outputAlphaFile.is_open()) {
		outputAlphaFile.precision(outputPrecision);
		outputAlphaFile << bestCLL << "\t" << bestICL << endl;
		engine.printDoubleVector2TextFile(bestAlpha, outputAlphaFile);
		outputAlphaFile.close();
	} else
		cout << "Can not open the file " << outputAlphaFilename << endl;

	for (int t = 0; t < engine.getNumOfParameterSets(); t++) {
		stringstream sT;
		sT << t;
		string outputPiFilename = outputDirectory + "/PI_" + sT.str() + "_"
				+ postfix.str() + ".txt";
		engine.printDoubleMatrix2TextFile(bestPIs(t), outputPiFilename.c_str(),
				outputPrecision);
	}

	stringstream bestTauPostfix;
	bestTauPostfix << engine.getMinNumOfClusters() << "_";
	bestTauPostfix << engine.getMaxNumOfEMIterations() << "_";
	bestTauPostfix << engine.getMaxNumOfESteps();

	string bestICLFilename = outputDirectory + "/bestICL_"
			+ bestTauPostfix.str() + ".txt";
	ifstream readBestICLFile;
	readBestICLFile.open(bestICLFilename.c_str());
	if (readBestICLFile.is_open()) {
		double currentICL;
		readBestICLFile >> currentICL;
		readBestICLFile.close();
		if (currentICL < bestICL) {
			string ouputTauFilename = outputDirectory + "/bestTau_"
					+ bestTauPostfix.str() + ".txt";
			engine.printDoubleMatrix2TextFile(bestTau, ouputTauFilename.c_str(),
					outputPrecision);
			// overwrite old ICL
			ofstream writeBestICLFile;
			writeBestICLFile.open(bestICLFilename.c_str());
			if (writeBestICLFile.is_open()) {
				writeBestICLFile.precision(outputPrecision);
				writeBestICLFile << bestICL;
				writeBestICLFile.close();
			} else
				cout << "Can not open the file " << bestICLFilename << endl;
		}
	} else {
		ofstream writeBestICLFile;
		writeBestICLFile.open(bestICLFilename.c_str());
		if (writeBestICLFile.is_open()) {
			writeBestICLFile.precision(outputPrecision);
			writeBestICLFile << bestICL;
			writeBestICLFile.close();
			string ouputTauFilename = outputDirectory + "/bestTau_"
					+ bestTauPostfix.str() + ".txt";
			engine.printDoubleMatrix2TextFile(bestTau, ouputTauFilename.c_str(),
					outputPrecision);
		} else
			cout << "Can not open the file " << bestICLFilename << endl;
	}
}

void estimateParametersWithConvergenceDetails(const char* confFile,
		const char* outputDir, unsigned int numOfClusters, unsigned int fromRun,
		unsigned int toRun) {

	cout << "Testing Signed Reciprocity Model " << endl;
	cout << "confFile is " << confFile << endl;
	cout << "outputDirectory is " << outputDir << endl;
	cout << "numOfClusters is " << numOfClusters << endl;

	EMEngine engine = EMEngine();
	engine.configure(confFile, "\t");
	engine.setMinNumOfClusters(numOfClusters);
	engine.setMaxNumOfClusters(numOfClusters);
	engine.setNumOfEMRuns(toRun - fromRun + 1);
	set_seed(fromRun, toRun);
	cout << "Finished configuration" << endl;

	DoubleMatrix bestTau = DoubleMatrix(engine.getNumOfVertices(),
			engine.getMinNumOfClusters());
	DoubleVector bestAlpha = DoubleVector(engine.getMinNumOfClusters());
	VectorOfDoubleMatrices bestPIs = VectorOfDoubleMatrices(
			engine.getNumOfParameterSets());
	for (int t = 0; t < engine.getNumOfParameterSets(); t++) {
		bestPIs(t) = DoubleMatrix(engine.getMinNumOfClusters(),
				engine.getMinNumOfClusters());
	}
	double bestCLL = 0;
	double bestICL = 0;

	engine.getBestEMRunConvergenceDetails(engine.getMinNumOfClusters(),
			engine.getMaxNumOfEMIterations(), engine.getMaxNumOfESteps(),
			engine.getTauPrecision(), engine.getParamPrecision(),
			engine.getNumOfEMRuns(), bestTau, bestAlpha, bestPIs, bestCLL,
			bestICL);

	string outputDirectory = string(outputDir);
	unsigned int outputPrecision = 10;

	stringstream postfix;
	postfix << engine.getMinNumOfClusters() << "_";
	postfix << fromRun << "_";
	postfix << toRun << "_";
	postfix << engine.getMaxNumOfEMIterations() << "_";
	postfix << engine.getMaxNumOfESteps();

	string outputAlphaFilename = outputDirectory + "/alpha_" + postfix.str()
			+ ".txt";
	ofstream outputAlphaFile;
	outputAlphaFile.open(outputAlphaFilename.c_str());
	if (outputAlphaFile.is_open()) {
		outputAlphaFile.precision(outputPrecision);
		outputAlphaFile << bestCLL << "\t" << bestICL << endl;
		engine.printDoubleVector2TextFile(bestAlpha, outputAlphaFile);
		outputAlphaFile.close();
	} else
		cout << "Can not open the file " << outputAlphaFilename << endl;

	for (int t = 0; t < engine.getNumOfParameterSets(); t++) {
		stringstream sT;
		sT << t;
		string outputPiFilename = outputDirectory + "/PI_" + sT.str() + "_"
				+ postfix.str() + ".txt";
		engine.printDoubleMatrix2TextFile(bestPIs(t), outputPiFilename.c_str(),
				outputPrecision);
	}

	stringstream bestTauPostfix;
	bestTauPostfix << engine.getMinNumOfClusters() << "_";
	bestTauPostfix << engine.getMaxNumOfEMIterations() << "_";
	bestTauPostfix << engine.getMaxNumOfESteps();

	string bestICLFilename = outputDirectory + "/bestICL_"
			+ bestTauPostfix.str() + ".txt";
	ifstream readBestICLFile;
	readBestICLFile.open(bestICLFilename.c_str());
	if (readBestICLFile.is_open()) {
		double currentICL;
		readBestICLFile >> currentICL;
		readBestICLFile.close();
		if (currentICL < bestICL) {
			string ouputTauFilename = outputDirectory + "/bestTau_"
					+ bestTauPostfix.str() + ".txt";
			engine.printDoubleMatrix2TextFile(bestTau, ouputTauFilename.c_str(),
					outputPrecision);
			// overwrite old ICL
			ofstream writeBestICLFile;
			writeBestICLFile.open(bestICLFilename.c_str());
			if (writeBestICLFile.is_open()) {
				writeBestICLFile.precision(outputPrecision);
				writeBestICLFile << bestICL;
				writeBestICLFile.close();
			} else
				cout << "Can not open the file " << bestICLFilename << endl;
		}
	} else {
		ofstream writeBestICLFile;
		writeBestICLFile.open(bestICLFilename.c_str());
		if (writeBestICLFile.is_open()) {
			writeBestICLFile.precision(outputPrecision);
			writeBestICLFile << bestICL;
			writeBestICLFile.close();
			string ouputTauFilename = outputDirectory + "/bestTau_"
					+ bestTauPostfix.str() + ".txt";
			engine.printDoubleMatrix2TextFile(bestTau, ouputTauFilename.c_str(),
					outputPrecision);
		} else
			cout << "Can not open the file " << bestICLFilename << endl;
	}
}

void bootstrapSignedNetworksForStdErrors(const char* confFile,
		const char* inputDir, const char* outputDir, int numOfVertices,
		int numOfClusters, int fromRun, int toRun, int numOfEMIterations,
		int numOfESteps) {

	ParametricBoostrapSampler* sampler =
			new EfficientSignedReciprocityParametricBoostrapSampler();

	unsigned int outputPrecision = 10;
	DoubleVector bestAlpha = DoubleVector(numOfClusters);
	VectorOfDoubleMatrices bestPIs = VectorOfDoubleMatrices(6);
	for (int t = 0; t < 6; t++) {
		bestPIs(t) = DoubleMatrix(numOfClusters, numOfClusters);
	}

	cout << "Working on setting parameters!" << endl;
	stringstream postfix;
	postfix << numOfClusters << "_";
	postfix << fromRun << "_";
	postfix << toRun << "_";
	postfix << numOfEMIterations << "_";
	postfix << numOfESteps;

	string inputDirectory = string(inputDir);
	string inputAlphaFilename = inputDirectory + "/alpha_" + postfix.str()
			+ ".txt";
	cout << "Read from the alpha file " << inputAlphaFilename << endl;
	sampler->readDoubleVectorFromTextFileWithFirstLineRemove(bestAlpha,
			inputAlphaFilename.c_str());
	for (int t = 0; t < 6; t++) {
		stringstream sT;
		sT << t;
		string inputPiFilename = inputDirectory + "/PI_" + sT.str() + "_"
				+ postfix.str() + ".txt";
		cout << "Read from the PI file " << inputPiFilename << endl;
		sampler->readDoubleMatrixFromTextFile((DoubleMatrix&) bestPIs(t),
				inputPiFilename.c_str());
	}
	printAlpha(bestAlpha);
	printSignedPIs(bestPIs);
	sampler->setParameters(numOfVertices, numOfClusters, bestAlpha, bestPIs);
	cout << "Done with setting parameters!" << endl;

	cout << "Working on setting configuration!" << endl;
	string confDirectory = string(confFile);
	sampler->configure(confDirectory.c_str(), "\t");
	set_seed(sampler->getSeed1(), sampler->getSeed2());
	cout << "Done with setting configuration!" << endl;

	cout << "Working on bootstrapping samples!" << endl;
	string outputDirectory = string(outputDir);
	sampler->obtainThetaBootstrapSamples(outputDirectory.c_str(),
			outputPrecision);
	cout << "Done with bootstrapping samples!" << endl;

	if (sampler != NULL)
		delete sampler;
}

void bootstrapBinaryNetworksForStdErrors(const char* confFile,
		const char* inputDir, const char* outputDir, int numOfVertices,
		int numOfClusters, int fromRun, int toRun, int numOfEMIterations,
		int numOfESteps) {

	ParametricBoostrapSampler* sampler =
			new EfficientBinaryReciprocityParametricBoostrapSampler();

	unsigned int outputPrecision = 10;
	DoubleVector bestAlpha = DoubleVector(numOfClusters);
	VectorOfDoubleMatrices bestPIs = VectorOfDoubleMatrices(3);
	for (int t = 0; t < 3; t++) {
		bestPIs(t) = DoubleMatrix(numOfClusters, numOfClusters);
	}

	cout << "Working on setting parameters!" << endl;
	stringstream postfix;
	postfix << numOfClusters << "_";
	postfix << fromRun << "_";
	postfix << toRun << "_";
	postfix << numOfEMIterations << "_";
	postfix << numOfESteps;

	string inputDirectory = string(inputDir);
	string inputAlphaFilename = inputDirectory + "/alpha_" + postfix.str()
			+ ".txt";
	cout << "Reading from the alpha file " << inputAlphaFilename << endl;
	sampler->readDoubleVectorFromTextFileWithFirstLineRemove(bestAlpha,
			inputAlphaFilename.c_str());
	for (int t = 0; t < 3; t++) {
		stringstream sT;
		sT << t;
		string inputPiFilename = inputDirectory + "/PI_" + sT.str() + "_"
				+ postfix.str() + ".txt";
		cout << "Reading from the PI file " << inputPiFilename << endl;
		sampler->readDoubleMatrixFromTextFile((DoubleMatrix&) bestPIs(t),
				inputPiFilename.c_str());
	}
	printAlpha(bestAlpha);
	printBinaryPIs(bestPIs);
	sampler->setParameters(numOfVertices, numOfClusters, bestAlpha, bestPIs);
	cout << "Done with setting parameters!" << endl;

	cout << "Working on setting configuration!" << endl;
	string confDirectory = string(confFile);
	sampler->configure(confDirectory.c_str(), "\t");
	set_seed(sampler->getSeed1(), sampler->getSeed2());
	cout << "Done with setting configuration!" << endl;

	cout << "Working on bootstrapping samples!" << endl;
	string outputDirectory = string(outputDir);
	sampler->obtainThetaBootstrapSamples(outputDirectory.c_str(),
			outputPrecision);
	cout << "Done with bootstrapping samples!" << endl;

	if (sampler != NULL)
		delete sampler;
}

int main(int argc, char *argv[]) {

	switch (atoi(argv[1])) {
	case 0:
		// const char* confFile: the path to the configuration file
		// const char* outputDir: the path to the directory where the output should be
		// unsigned int numOfClusters: the number of clusters
		// unsigned int fromRun: the beginning run number which is used for setting random seed 1
		// unsigned int toRun: the ending run number which is used for setting random seed 2
		// The last two parameters allow for paralleling the runs on LION clusters.
		// For example, we can specify 1 to 10, and 11 to 20 to cut 20 runs into 2 parallel sections.
		estimateParameters(argv[2], argv[3], atoi(argv[4]), atoi(argv[5]),
				atoi(argv[6]));
		break;
	case 1:
		// const char* confFile: the path to the configuration file
		// const char* outputDir: the path to the directory where the output should be
		// unsigned int numOfClusters: the number of clusters
		// unsigned int fromRun: the beginning run number which is used for setting random seed 1
		// unsigned int toRun: the ending run number which is used for setting random seed 2
		// The last two parameters allow for paralleling the runs on LION clusters.
		// For example, we can specify 1 to 10, and 11 to 20 to cut 20 runs into 2 parallel sections.
		estimateParametersWithConvergenceDetails(argv[2], argv[3],
				atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));
		break;
	case 2:
		//	const char* confFile: the path to the configuration file
		//	const char* inputDir: the path to the directory where model parameters are placed
		//	const char* outputDir: the path to the output directory
		//	int numOfVertices: the number of nodes
		//	int numOfClusters: the number of clusters
		//	int fromRun: the beginning run number which was set to generate model parameters based on data
		//	int toRun: the ending run number which was set to generate model parameters based on data
		//	int numOfEMIterations: the number of EM iterations which was set to generate model parameters based on data
		//	int numOfESteps: the number of E steps which was set to generate model parameters based on data
		//  The last for parameters should be define so that the program can look for the correct model parameter files
		//  which are supposed to produce by the first option.
		bootstrapSignedNetworksForStdErrors(argv[2], argv[3], argv[4],
				atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8]),
				atoi(argv[9]), atoi(argv[10]));
		break;
	case 3:
		//	const char* confFile: the path to the configuration file
		//	const char* inputDir: the path to the directory where model parameters are placed
		//	const char* outputDir: the path to the output directory
		//	int numOfVertices: the number of nodes
		//	int numOfClusters: the number of clusters
		//	int fromRun: the beginning run number which was set to generate model parameters based on data
		//	int toRun: the ending run number which was set to generate model parameters based on data
		//	int numOfEMIterations: the number of EM iterations which was set to generate model parameters based on data
		//	int numOfESteps: the number of E steps which was set to generate model parameters based on data
		//  The last for parameters should be define so that the program can look for the correct model parameter files
		//  which are supposed to produce by the first option.
		bootstrapBinaryNetworksForStdErrors(argv[2], argv[3], argv[4],
				atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8]),
				atoi(argv[9]), atoi(argv[10]));
		break;
	default:
		cout << "Please enter another first parameter!!!" << endl;
		cout
				<< "0: Run the EM method to estimate parameters WITHOUT convergence details"
				<< endl;
		cout
				<< "1: Run the EM method to estimate parameters WITH convergence details"
				<< endl;
		cout << "2: Run the bootstrap method to estimate standard errors for "
				<< "signed networks" << endl;
		cout << "3: Run the bootstrap method to estimate standard errors for "
				<< "binary networks" << endl;
		break;
	}

	return 0;
}

