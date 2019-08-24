/*
 * utils.h
 *
 *  Created on: Mar 27, 2011
 *      Author: duyvu
 */

#ifndef UTILS_H_
#define UTILS_H_

#include "model/DataTypes.h"

#include <iostream>

using namespace std;
using namespace mixrerg;

void printDoubleMatrix(const string& message, DoubleMatrix& matrix) {
	cout << message << endl;
	for (unsigned int i = 0; i < matrix.size1(); i++) {
		for (unsigned int j = 0; j < matrix.size2(); j++)
			cout << matrix(i, j) << "  ";
		cout << endl;
	}
}

void printDoubleVector(const string& message, DoubleVector& vector) {
	cout << message << endl;
	for (unsigned int i = 0; i < vector.size(); i++)
		cout << vector(i) << endl;
}

void printDoubleMatrix(const char* message, DoubleMatrix& matrix) {
	cout << message << endl;
	for (unsigned int i = 0; i < matrix.size1(); i++) {
		for (unsigned int j = 0; j < matrix.size2(); j++)
			cout << matrix(i, j) << "  ";
		cout << endl;
	}
}

void printDoubleVector(const char* message, DoubleVector& vector) {
	cout << message << endl;
	for (unsigned int i = 0; i < vector.size(); i++)
		cout << vector(i) << endl;
}

void printAlpha(const DoubleVector& alpha) {
	cout.precision(6);
	cout << "Printing ALPHA: " << endl;
	cout << alpha(0);
	unsigned int numOfClasses = alpha.size();
	for (unsigned int k = 1; k < numOfClasses; k++)
		cout << ",  " << alpha(k);
	cout << endl;
}

void printBinaryPIs(const VectorOfDoubleMatrices& PIs) {
	cout.precision(6);
	unsigned int numOfClasses = PIs(0).size1();
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

}

void printSignedPIs(const VectorOfDoubleMatrices& PIs) {
	cout.precision(6);
	unsigned int numOfClasses = PIs(0).size1();
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

#endif /* UTILS_H_ */
