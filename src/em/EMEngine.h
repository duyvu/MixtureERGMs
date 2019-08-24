/*
 * EMEngine.h
 *
 *  Created on: Jan 19, 2010
 *      Author: duyvu
 */

#ifndef EMENGINE_H_
#define EMENGINE_H_

#include "model/ReciprocityModel.h"

#include "util/StringTokenizer.h"

#include "iostream"
#include "fstream"

using namespace std;

namespace mixrerg {

class EMEngine {
protected:
	ReciprocityModel* model;
	int numOfEMRuns;
	int maxNumOfEMIterations;
	int maxNumOfESteps;
	double minTau;
	double minAlpha;
	double minPi;
	double tauPrecision;
	double paramPrecision;
	int minNumOfClusters;
	int maxNumOfClusters;
	// the user needs to retrieve these values and set up the state of random generator.
	int seed1;
	int seed2;
public:
	EMEngine();
	virtual ~EMEngine();

	virtual void configure(const char* configurationFile, const char* delim);
	void setData(const char* dataFile, int fileFormat, const char* delim);
	void setData(const NetworkSparseMatrix& networkData);

	void runEM(const int & _maxNumOfEMIterations, const int & _maxNumOfESteps,
			const double & _tauPrecision, const double & _paramPrecision);
	void runEMConvergenceDetails(const int & _maxNumOfEMIterations,
			const int & _maxNumOfESteps, const double & _tauPrecision,
			const double & _paramPrecision);
	void runEM();

	void runInitializedEM(int numOfClusters, const char* tauFile,
			DoubleMatrix& bestTau, DoubleVector& bestAlpha,
			VectorOfDoubleMatrices& bestPIs);

	void getBestEMRun(const int _numOfClusters,
			const int & _maxNumOfEMIterations, const int & _maxNumOfESteps,
			const double & _tauPrecision, const double & _paramPrecision,
			const int & _numOfEMRuns, DoubleMatrix& bestTau,
			DoubleVector& bestAlpha, VectorOfDoubleMatrices& bestPIs,
			double& bestCLL, double& bestICL);
	void getBestEMRunConvergenceDetails(const int _numOfClusters,
			const int & _maxNumOfEMIterations, const int & _maxNumOfESteps,
			const double & _tauPrecision, const double & _paramPrecision,
			const int & _numOfEMRuns, DoubleMatrix& bestTau,
			DoubleVector& bestAlpha, VectorOfDoubleMatrices& bestPIs,
			double& bestCLL, double& bestICL);
	void getBestEMRun(const int _numOfClusters, DoubleMatrix& bestTau,
			DoubleVector& bestAlpha, VectorOfDoubleMatrices& bestPIs,
			double& bestCLL, double& bestICL);

	double getParamPrecision() const {
		return paramPrecision;
	}

	double getTauPrecision() const {
		return tauPrecision;
	}

	void setParamPrecision(double paramPrecision) {
		this->paramPrecision = paramPrecision;
	}

	void setTauPrecision(double tauPrecision) {
		this->tauPrecision = tauPrecision;
	}

	int getNumOfVertices() const {
		if (model != NULL)
			return model->getNumOfVertices();
		else
			return -1;
	}

	int getNumOfEdges() const {
		if (model != NULL)
			return model->getNumOfEdges();
		else
			return -1;
	}

	int getNumOfClasses() const {
		if (model != NULL)
			return model->getNumOfClasses();
		else
			return -1;
	}

	int getNumOfParameterSets() const {
		if (model != NULL)
			return model->getNumOfParameterSets();
		else
			return -1;
	}

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

	int getNumOfEMRuns() const {
		return numOfEMRuns;
	}

	void setNumOfEMRuns(int numOfEMRuns) {
		this->numOfEMRuns = numOfEMRuns;
	}

	int getMaxNumOfEMIterations() const {
		return maxNumOfEMIterations;
	}

	void setMaxNumOfEMIterations(int maxNumOfEMIterations) {
		this->maxNumOfEMIterations = maxNumOfEMIterations;
	}

	int getMaxNumOfESteps() const {
		return maxNumOfESteps;
	}

	void setMaxNumOfESteps(int maxNumOfESteps) {
		this->maxNumOfESteps = maxNumOfESteps;
	}

	double getMinTau() const {
		return minTau;
	}

	void setMinTau(double minTau) {
		this->minTau = minTau;
	}

	double getMinAlpha() const {
		return minAlpha;
	}

	void setMinAlpha(double minAlpha) {
		this->minAlpha = minAlpha;
	}

	double getMinPi() const {
		return minPi;
	}

	void setMinPi(double minPi) {
		this->minPi = minPi;
	}

	int getMinNumOfClusters() const {
		return minNumOfClusters;
	}

	void setMinNumOfClusters(int minNumOfClusters) {
		this->minNumOfClusters = minNumOfClusters;
	}

	int getMaxNumOfClusters() const {
		return maxNumOfClusters;
	}

	void setMaxNumOfClusters(int maxNumOfClusters) {
		this->maxNumOfClusters = maxNumOfClusters;
	}

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
};

}

#endif /* EMENGINE_H_ */
