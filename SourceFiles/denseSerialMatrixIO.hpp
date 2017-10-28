#ifndef denseSerialMatrixIO_HPP
#define denseSerialMatrixIO_HPP

#include "PlayaSerialVectorType.hpp"
#include "PlayaDenseSerialMatrix.hpp"
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#include "SundanceExpr.hpp"
#include "PlayaSerialEpetraAdapter.hpp" // Allows me to convert between Epetra and Serial Vectors
#include "SundanceDiscreteSpace.hpp"
#include "SundanceDiscreteFunction.hpp"


using std::vector;
using std::ofstream;
using std::ifstream;
using namespace Playa;
using namespace PlayaExprTemplates;
using Sundance::Expr;
using Sundance::DiscreteSpace;

/** 
 * This function writes the elements of the DenseSerialMatrix pointed to by A 
 * to filename row-wise
 */
void writeDenseSerialMatrix(const RCP<DenseSerialMatrix>& A, string filename);

/** 
 * This function reads the data in filename as the entries into a DenseSerialMatrix row-wise 
 */
void readDenseSerialMatrix(RCP<DenseSerialMatrix>& A, string filename);

/**
 * This function reads the data in filename as the entries into a square matrix stored 
 * as an Expr row-wise
 */
void readListExprMatrix(Expr &A, int N, string filename);

/**
 * This function returns a Vector<double> that is the mean of the columns of the matrix A
 */
Vector<double> meanVector(const LinearOperator<double> A);

/**
 * This function assumes that each column of the DenseSerialMatrix reprents velocity 
 * and/or pressure data for a fluid in a spatial domain. It will return a spatial-only
 * dependent function which is calculated as the time mean over all nodes of the mesh
 */
Expr timeMeanFunctionGenerator(const LinearOperator<double> A, const DiscreteSpace ds);

/**
 * This function returns a matrix for which the mean of the columns has been
 * subtracted from each column for the supplied matrix
 */
LinearOperator<double> generateFluctuationMatrix(const LinearOperator<double> A);


/**
 * This function takes in an array of pointers to DenseSerialMatrix objects (all with the
 * same number of rows, not necessarily the same number of columns) and compiles them
 * into a single LinearOperator<double>
 */
LinearOperator<double> matrixAssembly(const Array<RCP<DenseSerialMatrix>> library);


#endif
