/*
 * DataTypes.h
 *
 *  Created on: Jan 17, 2010
 *      Author: duyvu
 */

#ifndef DATATYPES_H_
#define DATATYPES_H_

#include <cmath>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector_of_vector.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>

using namespace boost::numeric::ublas;

namespace mixrerg {

#define BINARY_GRAPH 0
#define MM_BINARY_GRAPH 1

#define SIGNED_GRAPH 2
#define MM_SIGNED_GRAPH 3

#define BOOTSTRAP_BINARY_GRAPH 4
#define BOOTSTRAP_SIGNED_GRAPH 5

#define SPARSE_FORMAT 0
#define MATRIX_FORMAT 1

#define BLANK_SPACE 0
#define TAB 1

#define MIN_TAU 1e-10
#define MIN_ALPHA 1e-6
#define MIN_PI 1e-6

#define TAU_PRECISION 1e-10
#define PI_PRECISION 1e-10

#define REL_LLH 1e-10

typedef compressed_matrix<signed short, row_major> NetworkSparseMatrix;
typedef generalized_vector_of_vector<signed int, row_major, vector<
		compressed_vector<signed int> > > SparseGVOCoordinateV;

typedef vector<unsigned int> UnsignedIntVector;
typedef vector<long> LongVector;
typedef vector<long double> DoubleVector;
typedef matrix<long double> DoubleMatrix;
typedef vector<DoubleMatrix> VectorOfDoubleMatrices;

typedef matrix<DoubleVector> MatrixOfDoubleVectors;
typedef vector<MatrixOfDoubleVectors> VectorOfMatrixOfDoubleVectors;

typedef matrix<double*> DoublePointerMatrix;

}

#endif /* DATATYPES_H_ */
