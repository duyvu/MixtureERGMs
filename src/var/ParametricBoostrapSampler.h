/*
 * ParametricBoostrapSampler.h
 *
 *  Created on: Oct 4, 2010
 *      Author: duyvu
 */

#ifndef PARAMETRICBOOSTRAPSAMPLER_H_
#define PARAMETRICBOOSTRAPSAMPLER_H_

#include "model/DataTypes.h"
#include "model/ReciprocityModel.h"
#include "util/StringTokenizer.h"

#include "iostream"
#include "fstream"

using namespace std;

namespace mixrerg {

class ParametricBoostrapSampler {
protected:
	unsigned int numOfVertices;
	unsigned int numOfClasses;
	unsigned int numOfPiParameterSets;
	DoubleVector alpha;
	VectorOfDoubleMatrices PIs;
	DoublePointerMatrix dyadProbs;

	// for the estimation task
	ReciprocityModel* model;
	int modelType;

	// for bootstrap
	unsigned int bootstrapSampleSize;

	// for the model precision
	double minTau;
	double minAlpha;
	double minPi;

	// for running the EM algorithm
	int numOfEMRuns;
	int maxNumOfEMIterations;
	int maxNumOfESteps;
	double tauPrecision;
	double paramPrecision;

	int seed1;
	int seed2;

	void runEM(const UnsignedIntVector& Z);
	void getBestEMRun(const UnsignedIntVector& Z, DoubleMatrix& bestTau,
			DoubleVector& bestAlpha, VectorOfDoubleMatrices& bestPIs,
			double& bestCLL, double& bestICL);
public:
	ParametricBoostrapSampler();
	virtual ~ParametricBoostrapSampler();

	// get the number of parameters theta
	virtual unsigned int getNumOfThetas() = 0;

	// for sampling networks
	virtual void setParameters(unsigned int _numOfVertices,
			unsigned int _numOfClasses, const DoubleVector& _alpha,
			const VectorOfDoubleMatrices& _PIs) = 0;
	virtual void normalizeDyadProbs();
	virtual void
	generateSample(NetworkSparseMatrix& Y, UnsignedIntVector& Z) = 0;
	void generateSamples(unsigned int numOfSamples,
			const string& outputDirectory, unsigned outputPrecision);
	void generateSamples(unsigned fromSampleIndex, unsigned toSampleIndex,
			const string& outputDirectory, unsigned outputPrecision);

	// for estimating parameters and get the bootstrap variances.
	void configure(const char* configurationFile, const char* delim);
	virtual void transformParameters(const VectorOfDoubleMatrices& PIs,
			VectorOfDoubleMatrices& Thetas) = 0;
	virtual unsigned int getNumberOfClusteringCoefficients() = 0;
	virtual void computeOverallClusteringCoefficients(
			DoubleVector& overallClusteringCoefficientSample,
			const DoubleVector& alpha, const VectorOfDoubleMatrices& PIs) = 0;
	void estimateSample(unsigned sampleIndex, const string& inputDirectory,
			const string& outputDirectory, unsigned int outputPrecision);
	void estimateSamples(unsigned fromSampleIndex, unsigned toSampleIndex,
			const string& inputDirectory, const string& outputDirectory,
			unsigned int outputPrecision);

	// summary the bootstrap samples
	void summarizeThetaSamples(unsigned fromSampleIndex,
			unsigned toSampleIndex, const string& inputDirectory,
			const string& outputDirectory, unsigned int outputPrecision);

	// integrate all tasks so that no intermediate output is generated
	void obtainThetaBootstrapSamples(const string& outputDirectory,
			const unsigned int outputPrecision);

	void printDoubleMatrix2TextFile(const DoubleMatrix& matrix,
			const char* filename, unsigned int outputPrecision) {
		ofstream outputFile;
		outputFile.open(filename);
		if (outputFile.is_open()) {
			outputFile.precision(outputPrecision);
			for (unsigned int k = 0; k < matrix.size1(); k++) {
				for (unsigned int l = 0; l < matrix.size2() - 1; l++)
					outputFile << matrix(k, l) << "\t";
				outputFile << matrix(k, matrix.size2() - 1);
				outputFile << endl;
			}
			outputFile.close();
		} else
			cout << "Can not open the file " << filename << endl;
	}

	void readDoubleMatrixFromTextFile(DoubleMatrix& matrix,
			const char* filename) {
		ifstream inputFile;
		inputFile.open(filename);
		matrix.clear();
		if (inputFile.is_open()) {
			string line;
			unsigned int rowIndex = 0;
			while (!inputFile.eof()) {
				getline(inputFile, line);
				unsigned int colIndex = 0;
				StringTokenizer tokenizer = StringTokenizer(line, "\t");
				while (tokenizer.hasMoreTokens()) {
					matrix(rowIndex, colIndex) = tokenizer.nextDoubleToken();
					colIndex++;
				}
				rowIndex++;
			}
			inputFile.close();
		} else
			cout << "Can not open the file " << filename << endl;
	}

	void printDoubleMatrix2TextFile(const DoubleMatrix& matrix,
			ofstream& outputFile) {
		for (unsigned int k = 0; k < matrix.size1(); k++) {
			for (unsigned int l = 0; l < matrix.size2() - 1; l++)
				outputFile << matrix(k, l) << "\t";
			outputFile << matrix(k, matrix.size2() - 1);
			outputFile << endl;
		}
	}

	void printDoubleVector2TextFile(const DoubleVector& vector,
			const char* filename, unsigned int outputPrecision) {
		ofstream outputFile;
		outputFile.open(filename);
		if (outputFile.is_open()) {
			outputFile.precision(outputPrecision);
			for (unsigned int k = 0; k < vector.size() - 1; k++)
				outputFile << vector(k) << "\t";
			outputFile << vector(vector.size() - 1);
			outputFile.close();
		} else
			cout << "Can not open the file " << filename << endl;
	}

	void printDoubleVector2TextFile(const DoubleVector& vector,
			ofstream& outputFile) {
		for (unsigned int k = 0; k < vector.size() - 1; k++) {
			outputFile << vector(k) << "\t";
		}
		outputFile << vector(vector.size() - 1);
	}

	void readDoubleVectorFromTextFileWithFirstLineRemove(DoubleVector& vector,
			const char* filename) {
		ifstream inputFile;
		inputFile.open(filename);
		vector.clear();
		if (inputFile.is_open()) {
			string line;
			getline(inputFile, line);
			//cout << "First line " << line << endl;
			getline(inputFile, line);
			//cout << "Second line " << line << endl;
			StringTokenizer tokenizer = StringTokenizer(line, "\t");
			unsigned colIndex = 0;
			while (tokenizer.hasMoreTokens()) {
				vector(colIndex++) = tokenizer.nextDoubleToken();
				//cout << colIndex - 1 << "\t" << vector(colIndex - 1) << endl;
			}
			inputFile.close();
		} else
			cout << "Can not open the file " << filename << endl;
	}

	void printAlpha(const DoubleVector& alpha);
	void printPIs(const VectorOfDoubleMatrices& PIs);
	int getSeed1() const {
		return seed1;
	}

	void setSeed1(int seed1) {
		this->seed1 = seed1;
	}

	int getSeed2() const {
		return seed2;
	}

	void setSeed2(int seed2) {
		this->seed2 = seed2;
	}

};

}

#endif /* PARAMETRICBOOSTRAPSAMPLER_H_ */
