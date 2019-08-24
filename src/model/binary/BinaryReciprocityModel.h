/*
 * BinaryReciprocityModel.h
 *
 *  Created on: Jan 17, 2010
 *      Author: duyvu
 */

#ifndef BINARYRECIPROCITYMODEL_H_
#define BINARYRECIPROCITYMODEL_H_

#include "model/ReciprocityModel.h"

namespace mixrerg {

class BinaryReciprocityModel: public mixrerg::ReciprocityModel {
protected:
	NetworkSparseMatrix stat10, stat01, stat11, stat00;

	DoubleMatrix tau, prevTau;
	DoubleVector alpha, prevAlpha;
	DoubleMatrix pi10, pi11, pi00; // pi11, pi00 are symmetric
	DoubleMatrix prevPi10, prevPi11;

	void updateTau(DoubleMatrix& new_tau, NetworkSparseMatrix& stat,
			const DoubleMatrix& tau, const DoubleMatrix& logPi,
			DoubleMatrix& temp);
	void updateTauByNegativeReflection(DoubleMatrix& new_tau,
			NetworkSparseMatrix& stat, const DoubleMatrix& tau,
			const DoubleMatrix& logPi, DoubleMatrix& temp);
	void updatePi(DoubleMatrix& pi, NetworkSparseMatrix& stat,
			const DoubleMatrix& tau, const DoubleMatrix& sumTaus);
public:
	BinaryReciprocityModel();
	virtual ~BinaryReciprocityModel();

	void setModel(int _numOfClasses, double _minTau, double _minAlpha,
			double _minPi);

	void resetParameters();
	void resetParameters(const UnsignedIntVector& Z);
	void resetParameters(const DoubleMatrix& initTau);
	void runFixedPointEstimationEStep();
	bool isTauSignificantlyChanged(double _tauPrecision);
	double getLargestTauChange();
	void runModelEstimationMStep();
	bool areAlphaPiSignificantlyChanged(double _paramPrecision);

	double getCompleteLogLikelihood();
	double getICL();
	void copyParameters(DoubleMatrix& bestTau, DoubleVector& bestAlpha,
			VectorOfDoubleMatrices& bestPIs);
	void copyPreviousParameters(DoubleMatrix& _prevTau,
			DoubleVector& _prevAlpha, VectorOfDoubleMatrices& _prevPIs);

	void printPIs();
	void printTau();
	void printAlpha();

};

}

#endif /* BINARYRECIPROCITYMODEL_H_ */
