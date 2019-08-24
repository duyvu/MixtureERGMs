/*
 * ReciprocityModel.h
 *
 *  Created on: Jan 17, 2010
 *      Author: duyvu
 */

#ifndef RECIPROCITYMODEL_H_
#define RECIPROCITYMODEL_H_

#include "DataTypes.h"

namespace mixrerg {

class ReciprocityModel {
protected:
	int numOfVertices;
	int numOfEdges;
	int numOfClasses;
	int numOfPiParameterSets;

	NetworkSparseMatrix networkMatrix;
	double minTau;
	double minAlpha;
	double minPi;

	double sumAbsoluteNetworkSparseMatrix(NetworkSparseMatrix& networkMatrix);
	void sumDoubleMatrixByRow(const DoubleMatrix& matrix, DoubleVector& vector);
	void sumDoubleMatrixByColumn(const DoubleMatrix& matrix,
			DoubleVector& vector);
	void logMatrix(DoubleMatrix& pi, DoubleMatrix& logPi);
	void logTransposedMatrix(DoubleMatrix& pi, DoubleMatrix& logPi);

	void normalizeLogTau2Tau(DoubleMatrix& tau, double minValue);
	void normalizeTau(DoubleMatrix& tau, double minValue);
	void normalizeVector(DoubleVector& vector, double minValue);

	/* Solve QP with one equality constraint and bounded constraints */
	void solveQP(const DoubleMatrix& m, const DoubleMatrix& s,
			DoubleMatrix& tau, double precision);
public:

	ReciprocityModel();
	virtual ~ReciprocityModel();

	void setData(const char* dataFile, int fileFormat, const char* delim);
	void setData(const NetworkSparseMatrix& networkData);

	virtual void setModel(int _numOfClasses, double _minTau, double _minAlpha,
			double _minPi) = 0;

	virtual void resetParameters() = 0;
	virtual void resetParameters(const UnsignedIntVector& Z) = 0;
	virtual void resetParameters(const DoubleMatrix& initTau) = 0;

	virtual void runFixedPointEstimationEStep() = 0;
	virtual bool isTauSignificantlyChanged(double _tauPrecision) = 0;
	virtual double getLargestTauChange() = 0;

	virtual void runModelEstimationMStep() = 0;
	virtual bool areAlphaPiSignificantlyChanged(double _paramPrecision) = 0;

	virtual double getCompleteLogLikelihood() = 0;
	virtual double getICL() = 0;
	virtual void copyParameters(DoubleMatrix& bestTau, DoubleVector& bestAlpha,
			VectorOfDoubleMatrices& bestPIs) = 0;
	virtual void copyPreviousParameters(DoubleMatrix& _prevTau,
			DoubleVector& _prevAlpha, VectorOfDoubleMatrices& _prevPIs) = 0;

	virtual void printPIs() = 0;
	virtual void printTau() = 0;
	virtual void printAlpha() = 0;

	void printDoubleMatrix(char* message, DoubleMatrix& matrix);
	void printDoubleVector(char* message, DoubleVector& vector);

	double getMinAlpha() const {
		return minAlpha;
	}

	double getMinPi() const {
		return minPi;
	}

	double getMinTau() const {
		return minTau;
	}

	void setMinAlpha(double minAlpha) {
		this->minAlpha = minAlpha;
	}

	void setMinPi(double minPi) {
		this->minPi = minPi;
	}

	void setMinTau(double minTau) {
		this->minTau = minTau;
	}

	int getNumOfVertices() const {
		return numOfVertices;
	}

	int getNumOfEdges() const {
		return numOfEdges;
	}

	int getNumOfClasses() const {
		return numOfClasses;
	}

	int getNumOfParameterSets() const {
		return numOfPiParameterSets;
	}

};

}

#endif /* RECIPROCITYMODEL_H_ */
