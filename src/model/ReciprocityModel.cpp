/*
 * ReciprocityModel.cpp
 *
 *  Created on: Jan 17, 2010
 *      Author: duyvu
 */

#include "ReciprocityModel.h"

#include "util/StringTokenizer.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <cfloat>

using namespace std;

namespace mixrerg {

ReciprocityModel::ReciprocityModel() {
	// TODO Auto-generated constructor stub
}

ReciprocityModel::~ReciprocityModel() {
	// TODO Auto-generated destructor stub
}

void ReciprocityModel::setData(const NetworkSparseMatrix& networkData) {
	networkMatrix = networkData;
	numOfVertices = networkMatrix.size1();
	numOfEdges = sumAbsoluteNetworkSparseMatrix(networkMatrix);
}

void ReciprocityModel::setData(const char* dataFile, int fileFormat,
		const char* delim) {

	// Entering data
	cout << "The network data file is " << dataFile << endl;

	string line;
	ifstream fileStream(dataFile);
	if (fileStream.is_open()) {

		// The first line is the number of vertices
		if (!fileStream.eof()) {
			getline(fileStream, line);
			StringTokenizer strtok = StringTokenizer(line, delim);
			numOfVertices = strtok.nextIntToken();
			cout << "The number of vertices is " << numOfVertices << endl;
			networkMatrix = NetworkSparseMatrix(numOfVertices, numOfVertices);

			// For sparse format matrix
			if (fileFormat == SPARSE_FORMAT) {
				cout << "SPARSE_FORMAT" << endl;
				while (!fileStream.eof()) {
					getline(fileStream, line);
					//cout << line << endl;
					StringTokenizer strtok = StringTokenizer(line, delim);
					//cout << strtok.countTokens() << endl;
					if (strtok.countTokens() == 3) {
						int i = strtok.nextIntToken();
						int j = strtok.nextIntToken();
						if (i == j) {
							cout << "Node " << i
									<< " has a self-edge. It will be removed."
									<< endl;
							continue;
						}
						signed short value = strtok.nextIntToken();
						//cout << "(" << i << ", " << j << ") = " << value << endl;
						if (networkMatrix(i, j) != 0)
							cout << "(" << i << ", " << j << ") = " << value
									<< " is duplicate!!!" << endl;
						else
							networkMatrix.insert_element(i, j, value);
					}
				}
			}

			// For normal matrix
			if (fileFormat == MATRIX_FORMAT) {
				cout << "MATRIX_FORMAT" << endl;
				int i = 0;
				while (!fileStream.eof()) {
					getline(fileStream, line);
					StringTokenizer strtok = StringTokenizer(line, delim);
					//cout << strtok.countTokens() << endl;
					if (strtok.countTokens() == numOfVertices) {
						for (int j = 0; j < numOfVertices; j++) {
							signed short value = strtok.nextIntToken();
							if (value != 0) {
								networkMatrix.insert_element(i, j, value);
							}
						}
					}
					i++;
				}
			}
			numOfEdges = sumAbsoluteNetworkSparseMatrix(networkMatrix);
			cout << "The number of edges is " << numOfEdges << endl;
		}

		fileStream.close();
	} else
		cout << "Unable to open file " << dataFile << endl;
}

double ReciprocityModel::sumAbsoluteNetworkSparseMatrix(
		NetworkSparseMatrix& netMatrix) {
	double result = 0;
	typedef NetworkSparseMatrix::iterator1 i1_t;
	typedef NetworkSparseMatrix::iterator2 i2_t;
	for (i1_t i1 = netMatrix.begin1(); i1 != netMatrix.end1(); ++i1) {
		for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
			result += abs(*i2);
	}
	return result;
}

void ReciprocityModel::normalizeTau(DoubleMatrix& tau, double minValue) {
	// normalize
	for (int i = 0; i < numOfVertices; i++) {
		double denominator = 0;
		for (int k = 0; k < numOfClasses; k++) {
			denominator += tau(i, k);
		}
		bool again = false;
		for (int k = 0; k < numOfClasses; k++) {
			tau(i, k) /= denominator;
			if (tau(i, k) < minValue) {
				tau(i, k) = minValue;
				again = true;
			}
		}
		if (again) {
			denominator = 0;
			for (int k = 0; k < numOfClasses; k++) {
				denominator += tau(i, k);
			}
			for (int k = 0; k < numOfClasses; k++)
				tau(i, k) /= denominator;
		}
	}
}

void ReciprocityModel::normalizeLogTau2Tau(DoubleMatrix& tau, double minValue) {
	int numOfVertices = tau.size1();
	int numOfClasses = tau.size2();

	const double logDoubleMax = log(DBL_MAX) - 1;

	for (int i = 0; i < numOfVertices; i++) {
		// First it holds the max value
		double slidingValue = tau(i, 0);
		for (int k = 1; k < numOfClasses; k++)
			if (slidingValue < tau(i, k))
				slidingValue = tau(i, k);
		// Now it actually holds the sliding value
		slidingValue = logDoubleMax - slidingValue;
		for (int k = 0; k < numOfClasses; k++)
			tau(i, k) += slidingValue;
	}

	// normalize
	for (int i = 0; i < numOfVertices; i++) {
		double denominator = 0;
		for (int k = 0; k < numOfClasses; k++) {
			tau(i, k) = exp(tau(i, k));
			denominator += tau(i, k);
		}
		bool again = false;
		for (int k = 0; k < numOfClasses; k++) {
			tau(i, k) /= denominator;
			if (tau(i, k) < minValue) {
				tau(i, k) = minValue;
				again = true;
			}
		}
		if (again) {
			denominator = 0;
			for (int k = 0; k < numOfClasses; k++) {
				denominator += tau(i, k);
			}
			for (int k = 0; k < numOfClasses; k++)
				tau(i, k) /= denominator;
		}
	}
}

void ReciprocityModel::normalizeVector(DoubleVector& vector, double minValue) {
	for (unsigned int k = 0; k < vector.size(); k++)
		if (vector(k) < minValue)
			vector(k) = minValue;
	double sumAlpha = 0;
	for (unsigned int k = 0; k < vector.size(); k++)
		sumAlpha += vector(k);
	for (unsigned int k = 0; k < vector.size(); k++)
		vector(k) /= sumAlpha;
}

void ReciprocityModel::logMatrix(DoubleMatrix& pi, DoubleMatrix& logPi) {
	for (unsigned int k = 0; k < pi.size1(); k++)
		for (unsigned int l = 0; l < pi.size2(); l++)
			logPi(k, l) = log(pi(k, l));
}

void ReciprocityModel::logTransposedMatrix(DoubleMatrix& pi,
		DoubleMatrix& logPi) {
	for (unsigned int k = 0; k < pi.size1(); k++)
		for (unsigned int l = 0; l < pi.size2(); l++)
			logPi(l, k) = log(pi(k, l));
}

void ReciprocityModel::sumDoubleMatrixByRow(const DoubleMatrix & matrix,
		DoubleVector & vector) {
	for (unsigned int j = 0; j < matrix.size2(); j++) {
		vector(j) = 0;
		for (unsigned int i = 0; i < matrix.size1(); i++)
			vector(j) += matrix(i, j);
	}
}

void ReciprocityModel::sumDoubleMatrixByColumn(const DoubleMatrix & matrix,
		DoubleVector & vector) {
	for (unsigned int i = 0; i < matrix.size1(); i++) {
		vector(i) = 0;
		for (unsigned int j = 0; j < matrix.size2(); j++)
			vector(i) += matrix(i, j);
	}
}

void ReciprocityModel::printDoubleMatrix(char* message, DoubleMatrix& matrix) {
	cout << message << endl;
	for (unsigned int i = 0; i < matrix.size1(); i++) {
		for (unsigned int j = 0; j < matrix.size2(); j++)
			cout << matrix(i, j) << "  ";
		cout << endl;
	}
}

void ReciprocityModel::printDoubleVector(char* message, DoubleVector& vector) {
	cout << message << endl;
	for (unsigned int i = 0; i < vector.size(); i++)
		cout << vector(i) << endl;
}

/* Solve QP with one equality constraint and bounded constraints */
void ReciprocityModel::solveQP(const DoubleMatrix& m, const DoubleMatrix& s,
		DoubleMatrix& tau, double precision) {

	if (m.size1() <= 0)
		return;

	unsigned int n = m.size2();

	bool J_k[n];
	bool J_lambda_a[n];

	bool J_lambda_k[n];
	bool J_lambda_k_a[n];
	bool J_lambda_k_b[n];

	double lambda_k;

	for (unsigned int rowIndex = 0; rowIndex < m.size1(); rowIndex++) {
		// Step 1:
		for (unsigned int j = 0; j < n; j++) {
			J_k[j] = true;
			J_lambda_a[j] = false;
		}

		bool done = false;
		// we never loop over n times.
		unsigned int count = 0;
		while (!done) {
			count++;

			// Step 2:
			double value1 = 0, value2 = 0;
			for (unsigned int j = 0; j < n; j++)
				if (J_k[j]) {
					value1 += 1 / m(rowIndex, j);
					value2 += s(rowIndex, j) / m(rowIndex, j);
				}
			lambda_k = (value2 - 2) / value1;

			// Step3:
			for (unsigned int j = 0; j < n; j++) {
				if (J_k[j]) {
					if (lambda_k >= s(rowIndex, j)) {
						J_lambda_k[j] = false;
						J_lambda_k_a[j] = true;
						J_lambda_k_b[j] = false;
					} else if (lambda_k
							> (-2 * m(rowIndex, j) + s(rowIndex, j))) {
						J_lambda_k[j] = true;
						J_lambda_k_a[j] = false;
						J_lambda_k_b[j] = false;
					} else {
						J_lambda_k[j] = false;
						J_lambda_k_a[j] = false;
						J_lambda_k_b[j] = true;
					}
				} else {
					J_lambda_k[j] = false;
					J_lambda_k_a[j] = false;
					J_lambda_k_b[j] = false;
				}
			}

			// Step 4:
			double delta = 0;
			double value3 = 0, value4 = 0;
			bool J_lambda_k_empty = true;
			for (unsigned int j = 0; j < n; j++) {
				delta += J_lambda_k_b[j];
				if (J_lambda_k[j]) {
					value3 += s(rowIndex, j) / m(rowIndex, j);
					value4 += 1 / m(rowIndex, j);
					J_lambda_k_empty = false;
				}
			}
			delta += (value3 / 2 - lambda_k * value4 / 2 - 1);

			if (fabs(delta) < precision || J_lambda_k_empty || count >= n) {
				for (unsigned int j = 0; j < n; j++) {
					if (J_lambda_a[j] || J_lambda_k_a[j])
						tau(rowIndex, j) = 0;
					else if (J_lambda_k_b[j])
						tau(rowIndex, j) = 1;
					else
						tau(rowIndex, j) = (s(rowIndex, j) - lambda_k) / (2
								* m(rowIndex, j));
				}
				done = true;
			} else if (delta > 0) {
				for (unsigned int j = 0; j < n; j++) {
					if (J_lambda_k_a[j]) {
						J_lambda_a[j] = true;
						J_k[j] = false;
					}
				}
			} else {
				for (unsigned int j = 0; j < n; j++) {
					if (J_lambda_k_b[j])
						tau(rowIndex, j) = 1;
					else
						tau(rowIndex, j) = 0;
				}
				done = true;
			}
		}
	}
}

}
