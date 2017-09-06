#ifndef denseSerialMatrixIO_HPP
#define denseSerialMatrixIO_HPP

#include "PlayaSerialVectorType.hpp"
#include "PlayaDenseSerialMatrix.hpp"
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#include "SundanceExpr.hpp"


using std::vector;
using std::ofstream;
using std::ifstream;
using namespace Playa;
using namespace PlayaExprTemplates;
using Sundance::Expr;

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


#endif
