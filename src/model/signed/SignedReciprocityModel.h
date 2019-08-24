/*
 * SignedReciprocityModel.h
 *
 *  Created on: Jan 17, 2010
 *      Author: duyvu
 */

#ifndef SIGNEDRECIPROCITYMODEL_H_
#define SIGNEDRECIPROCITYMODEL_H_

#include "model/ReciprocityModel.h"

namespace mixrerg {

class SignedReciprocityModel: public mixrerg::ReciprocityModel {
protected:

	NetworkSparseMatrix stat1a, stat1b, stat2a, stat2b, stat3a, stat3b, stat4m,
			stat4p, stat5;
	NetworkSparseMatrix statMinus, statPlus;

	DoubleMatrix tau, prevTau;
	DoubleVector alpha, prevAlpha;
	DoubleMatrix pi1, pi2, pi3, pi4m, pi4p, pi5;
	DoubleMatrix prevPi1, prevPi2, prevPi3, prevPi4m, prevPi4p;

	void updateTau(DoubleMatrix& new_tau, NetworkSparseMatrix& stat,
			const DoubleMatrix& tau, const DoubleMatrix& logPi,
			DoubleMatrix& temp);
	void updateTauByNegativeReflection(DoubleMatrix& new_tau,
			NetworkSparseMatrix& stat, const DoubleMatrix& tau,
			const DoubleMatrix& logPi, DoubleMatrix& temp);
	void updatePi(DoubleMatrix& pi, NetworkSparseMatrix& stat,
			const DoubleMatrix& tau, const DoubleMatrix& sumTaus);
public:
	SignedReciprocityModel();
	virtual ~SignedReciprocityModel();

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

#endif /* SIGNEDRECIPROCITYMODEL_H_ */
