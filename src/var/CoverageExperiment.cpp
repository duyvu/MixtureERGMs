/*
 * CoverageExperiment.cpp
 *
 *  Created on: Oct 13, 2010
 *      Author: duyvu
 */

#include "CoverageExperiment.h"

namespace mixrerg {

CoverageExperiment::CoverageExperiment() {
	// TODO Auto-generated constructor stub
	sampler = NULL;
}

CoverageExperiment::~CoverageExperiment() {
	// TODO Auto-generated destructor stub
	if (sampler == NULL)
		delete sampler;
}

void CoverageExperiment::getCoverageSample(unsigned int coverageSampleIndex) {

	cout << "coverageSampleIndex " << coverageSampleIndex << endl;

	NetworkSparseMatrix Y = NetworkSparseMatrix(numOfVertices, numOfVertices);
	UnsignedIntVector Z = UnsignedIntVector(numOfVertices);
	generateSample(Y, Z);

	DoubleMatrix bestTau = DoubleMatrix(numOfVertices, numOfClasses);
	DoubleVector bestAlpha = DoubleVector(numOfClasses);
	VectorOfDoubleMatrices bestPIs = VectorOfDoubleMatrices(
			numOfPiParameterSets);
	for (unsigned int t = 0; t < numOfPiParameterSets; t++)
		bestPIs(t) = DoubleMatrix(numOfClasses, numOfClasses);
	double bestCLL = 0;
	double bestICL = 0;

	engine.setData(Y);
	engine.getBestEMRun(numOfClasses, bestTau, bestAlpha, bestPIs, bestCLL,
			bestICL);

	sampler->setParameters(numOfVertices, numOfClasses, bestAlpha, bestPIs);
	stringstream sCoverageSampleIndex;
	sCoverageSampleIndex << outputDirectory << "/" << coverageSampleIndex;
	string newDirectory = "mkdir " + sCoverageSampleIndex.str();
	system(newDirectory.c_str());
	sampler->obtainThetaBootstrapSamples(sCoverageSampleIndex.str(),
			outputPrecision);
}

void CoverageExperiment::getCoverageSamples(
		unsigned int fromCoverageSampleIndex,
		unsigned int toCoverageSampleIndex) {

	NetworkSparseMatrix Y = NetworkSparseMatrix(numOfVertices, numOfVertices);
	UnsignedIntVector Z = UnsignedIntVector(numOfVertices);
	DoubleMatrix bestTau = DoubleMatrix(numOfVertices, numOfClasses);
	DoubleVector bestAlpha = DoubleVector(numOfClasses);
	VectorOfDoubleMatrices bestPIs = VectorOfDoubleMatrices(
			numOfPiParameterSets);
	for (unsigned int t = 0; t < numOfPiParameterSets; t++)
		bestPIs(t) = DoubleMatrix(numOfClasses, numOfClasses);
	double bestCLL = 0;
	double bestICL = 0;

	for (unsigned int coverageSampleIndex = fromCoverageSampleIndex; coverageSampleIndex
			<= toCoverageSampleIndex; coverageSampleIndex++) {

		cout << "coverageSampleIndex " << coverageSampleIndex << endl;

		generateSample(Y, Z);

		engine.setData(Y);
		engine.getBestEMRun(numOfClasses, bestTau, bestAlpha, bestPIs, bestCLL,
				bestICL);

		sampler->setParameters(numOfVertices, numOfClasses, bestAlpha, bestPIs);
		stringstream sCoverageSampleIndex;
		sCoverageSampleIndex << outputDirectory << "/" << coverageSampleIndex;
		string newDirectory = "mkdir " + sCoverageSampleIndex.str();
		system(newDirectory.c_str());
		sampler->obtainThetaBootstrapSamples(sCoverageSampleIndex.str(),
				outputPrecision);

	}
}

void CoverageExperiment::getCoverageSamples() {

	NetworkSparseMatrix Y = NetworkSparseMatrix(numOfVertices, numOfVertices);
	UnsignedIntVector Z = UnsignedIntVector(numOfVertices);
	DoubleMatrix bestTau = DoubleMatrix(numOfVertices, numOfClasses);
	DoubleVector bestAlpha = DoubleVector(numOfClasses);
	VectorOfDoubleMatrices bestPIs = VectorOfDoubleMatrices(
			numOfPiParameterSets);
	for (unsigned int t = 0; t < numOfPiParameterSets; t++)
		bestPIs(t) = DoubleMatrix(numOfClasses, numOfClasses);
	double bestCLL = 0;
	double bestICL = 0;

	for (unsigned int coverageSampleIndex = 0; coverageSampleIndex
			< coverageSampleSize; coverageSampleIndex++) {

		cout << "coverageSampleIndex " << coverageSampleIndex << endl;

		generateSample(Y, Z);

		engine.setData(Y);
		engine.getBestEMRun(numOfClasses, bestTau, bestAlpha, bestPIs, bestCLL,
				bestICL);

		sampler->setParameters(numOfVertices, numOfClasses, bestAlpha, bestPIs);
		stringstream sCoverageSampleIndex;
		sCoverageSampleIndex << outputDirectory << "/" << coverageSampleIndex;
		string newDirectory = "mkdir " + sCoverageSampleIndex.str();
		system(newDirectory.c_str());
		sampler->obtainThetaBootstrapSamples(sCoverageSampleIndex.str(),
				outputPrecision);
	}
}

void CoverageExperiment::getCoverageSampleUsingParameterFiles(
		unsigned int coverageSampleIndex) {

	DoubleVector bestAlpha = DoubleVector(numOfClasses);
	VectorOfDoubleMatrices bestPIs = VectorOfDoubleMatrices(
			numOfPiParameterSets);
	for (unsigned int t = 0; t < numOfPiParameterSets; t++)
		bestPIs(t) = DoubleMatrix(numOfClasses, numOfClasses);

	cout << "coverageSampleIndex " << coverageSampleIndex << endl;
	stringstream sCoverageSampleIndex;
	sCoverageSampleIndex << outputDirectory << "/" << coverageSampleIndex;

	string inputAlphaFilename = sCoverageSampleIndex.str() + "/alpha.txt";
	ifstream inputAlphaFile;
	inputAlphaFile.open(inputAlphaFilename.c_str());
	if (inputAlphaFile.is_open()) {
		string line;
		getline(inputAlphaFile, line);
		StringTokenizer tokenizer = StringTokenizer(line, "\t");
		unsigned int index = 0;
		while (tokenizer.hasMoreTokens()) {
			bestAlpha(index++) = tokenizer.nextDoubleToken();
		}
	} else
		cout << "Can not open the file " << inputAlphaFilename << endl;

	for (unsigned int indexOfPI = 0; indexOfPI < numOfPiParameterSets; indexOfPI++) {
		stringstream sIndexOfPI;
		sIndexOfPI << indexOfPI;
		string inputPIFilename = sCoverageSampleIndex.str() + "/PI"
				+ sIndexOfPI.str() + ".txt";
		sampler->readDoubleMatrixFromTextFile(bestPIs(indexOfPI),
				inputPIFilename.c_str());
	}

	sampler->setParameters(numOfVertices, numOfClasses, bestAlpha, bestPIs);
	sampler->obtainThetaBootstrapSamples(sCoverageSampleIndex.str(),
			outputPrecision);
}

void CoverageExperiment::getCoverageSamplesUsingParameterFiles(
		unsigned int fromCoverageSampleIndex,
		unsigned int toCoverageSampleIndex) {

	NetworkSparseMatrix Y = NetworkSparseMatrix(numOfVertices, numOfVertices);
	UnsignedIntVector Z = UnsignedIntVector(numOfVertices);
	DoubleMatrix bestTau = DoubleMatrix(numOfVertices, numOfClasses);
	DoubleVector bestAlpha = DoubleVector(numOfClasses);
	VectorOfDoubleMatrices bestPIs = VectorOfDoubleMatrices(
			numOfPiParameterSets);
	for (unsigned int t = 0; t < numOfPiParameterSets; t++)
		bestPIs(t) = DoubleMatrix(numOfClasses, numOfClasses);
	double bestCLL = 0;
	double bestICL = 0;

	for (unsigned int coverageSampleIndex = fromCoverageSampleIndex; coverageSampleIndex
			<= toCoverageSampleIndex; coverageSampleIndex++) {

		cout << "coverageSampleIndex " << coverageSampleIndex << endl;
		stringstream sCoverageSampleIndex;
		sCoverageSampleIndex << outputDirectory << "/" << coverageSampleIndex;

		string inputAlphaFilename = sCoverageSampleIndex.str() + "/alpha.txt";
		ifstream inputAlphaFile;
		inputAlphaFile.open(inputAlphaFilename.c_str());
		if (inputAlphaFile.is_open()) {
			string line;
			getline(inputAlphaFile, line);
			StringTokenizer tokenizer = StringTokenizer(line, "\t");
			unsigned int index = 0;
			while (tokenizer.hasMoreTokens()) {
				bestAlpha(index++) = tokenizer.nextDoubleToken();
			}
		} else
			cout << "Can not open the file " << inputAlphaFilename << endl;

		for (unsigned int indexOfPI = 0; indexOfPI < numOfPiParameterSets; indexOfPI++) {
			stringstream sIndexOfPI;
			sIndexOfPI << indexOfPI;
			string inputPIFilename = sCoverageSampleIndex.str() + "/PI"
					+ sIndexOfPI.str() + ".txt";
			sampler->readDoubleMatrixFromTextFile(bestPIs(indexOfPI),
					inputPIFilename.c_str());
		}
		sampler->setParameters(numOfVertices, numOfClasses, bestAlpha, bestPIs);
		sampler->obtainThetaBootstrapSamples(sCoverageSampleIndex.str(),
				outputPrecision);

	}
}

void CoverageExperiment::getCoverageSamplesUsingParameterFiles() {

	DoubleVector bestAlpha = DoubleVector(numOfClasses);
	VectorOfDoubleMatrices bestPIs = VectorOfDoubleMatrices(
			numOfPiParameterSets);
	for (unsigned int t = 0; t < numOfPiParameterSets; t++)
		bestPIs(t) = DoubleMatrix(numOfClasses, numOfClasses);

	for (unsigned int coverageSampleIndex = 0; coverageSampleIndex
			< coverageSampleSize; coverageSampleIndex++) {

		cout << "coverageSampleIndex " << coverageSampleIndex << endl;
		stringstream sCoverageSampleIndex;
		sCoverageSampleIndex << outputDirectory << "/" << coverageSampleIndex;

		string inputAlphaFilename = sCoverageSampleIndex.str() + "/alpha.txt";
		ifstream inputAlphaFile;
		inputAlphaFile.open(inputAlphaFilename.c_str());
		if (inputAlphaFile.is_open()) {
			string line;
			getline(inputAlphaFile, line);
			StringTokenizer tokenizer = StringTokenizer(line, "\t");
			unsigned int index = 0;
			while (tokenizer.hasMoreTokens()) {
				bestAlpha(index++) = tokenizer.nextDoubleToken();
			}
		} else
			cout << "Can not open the file " << inputAlphaFilename << endl;

		for (unsigned int indexOfPI = 0; indexOfPI < numOfPiParameterSets; indexOfPI++) {
			stringstream sIndexOfPI;
			sIndexOfPI << indexOfPI;
			string inputPIFilename = sCoverageSampleIndex.str() + "/PI"
					+ sIndexOfPI.str() + ".txt";
			sampler->readDoubleMatrixFromTextFile(bestPIs(indexOfPI),
					inputPIFilename.c_str());
		}

		sampler->setParameters(numOfVertices, numOfClasses, bestAlpha, bestPIs);
		sampler->obtainThetaBootstrapSamples(sCoverageSampleIndex.str(),
				outputPrecision);
	}
}

}
