/* @HEADER@ */
// ************************************************************************
// 
//                 Playa: Programmable Linear Algebra
//                 Copyright 2012 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
// 
// Author: Simon Rush
// Date: 10/04/16
/* @HEADER@ */

#include "PlayaDenseSerialMatrix.hpp"
#include "PlayaDenseMatrixMatrixProduct.hpp"
#include "PlayaSerialVector.hpp"
#include "PlayaExceptions.hpp"
//#include "EpetraExt_MatrixMatrix.h"
#include "Teuchos_LAPACK.hpp" // Needed for dgemm
#include "PlayaVectorType.hpp" // For VectorType, VectorSpace
#include "PlayaSerialVectorType.hpp" // For SerialVectorType

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#endif

/*! \cond */
extern "C"
{
void dgemm_(  char* TRANSA,
	      char* TRANSB,
	      int*   M,
	      int* 	 N,
	      int* 	 K,
              double*      ALPHA,
double* 	A,
	      int* 	LDA,
double* 	B,
	      int* 	LDB,
       double*          BETA,
double* 	C,
	      int* 	LDC);
}
/*! \endcond */

namespace Playa
{
using namespace Teuchos;

// This function performs the multiplication dA, where A matrix, and d is a diagonal matrix stored
// as a Playa::Vector. Since we are multiplying on the left, each element in d scales
// the corresponding row in A
LinearOperator<double> denseLeftScale(
  const Vector<double>& d,
  const LinearOperator<double>& A)
{
// Determine if A is a transpose
	const Playa::SimpleTransposedOp<double>* isAtrans = dynamic_cast<const Playa::SimpleTransposedOp<double>* >(A.ptr().get());

	// Extract the underlying PlayaDenseSerialMatrix. Type checking is done
	// by rcp_dynamic_cast, so we need no error checking here.
	RCP<Playa::DenseSerialMatrix> APtr;
	if(isAtrans)
		APtr = rcp_dynamic_cast<DenseSerialMatrix>(isAtrans->op().ptr());	
	else
		APtr = rcp_dynamic_cast<DenseSerialMatrix>(A.ptr());
	

	// Get a pointer to the DenseSerialMatrix's data_
	const double* A0_ptr = APtr->dataPtr();

	// Scale the return value
	RCP<DenseSerialMatrix> MPtr;
	LinearOperator<double> rtnMtx;
	double* rtn;
	if(!isAtrans) // A was given to denseLeftScale
	{
		TEUCHOS_TEST_FOR_EXCEPTION(APtr->numRows() != d.dim(), std::runtime_error, "The " << APtr->numRows() << "x" << APtr->numCols() << " matrix A was given with the vector d with dimension " << d.dim() );

		MPtr = rcp(new DenseSerialMatrix(A.domain(), A.range()));
/*
 * This is necessary so that we can return our matrix as a LinearOperator.
 * Both LinearOperator and DenseSerialMatrix inherit from LinearOperatorBase
 */
		rtnMtx = rcp_dynamic_cast<LinearOperatorBase<double> >(MPtr);
		rtn = MPtr->dataPtr();
		for(int row = 0; row<APtr->numRows(); row++)
		{
			for(int col = 0; col<APtr->numCols(); col++)
			{
				rtn[row+APtr->numRows()*col] = d[row]*A0_ptr[row+APtr->numRows()*col];
			}
		}
	}
	else // A.transpose() was given to denseLeftScale
	{
		TEUCHOS_TEST_FOR_EXCEPTION(APtr->numCols() != d.dim(), std::runtime_error, "The " << APtr->numCols() << "x" << APtr->numRows() << " matrix A^T was given with the vector d with dimension " << d.dim() );

		MPtr = rcp(new DenseSerialMatrix(isAtrans->op().range(), isAtrans->op().domain()));
/*
 * This is necessary so that we can return our matrix as a LinearOperator.
 * Both LinearOperator and DenseSerialMatrix inherit from LinearOperatorBase
 */
		rtnMtx = rcp_dynamic_cast<LinearOperatorBase<double> >(MPtr);
		rtn = MPtr->dataPtr();
		for(int row = 0; row<APtr->numCols(); row++)
		{
			for(int col = 0; col<APtr->numRows(); col++)
			{
				rtn[row+APtr->numCols()*col] = d[row]*A0_ptr[col+APtr->numRows()*row];
			}
		}
	}	

 	return rtnMtx; 
}

// This function performs the multiplication Ad, where A matrix, and d is a diagonal matrix stored
// as a Playa::Vector. Since we are multiplying on the right, each element in d scales
// the corresponding column in A
LinearOperator<double> denseRightScale(
  const LinearOperator<double>& A,
  const Vector<double>& d)
{


	// Determine if A is a transpose
	const Playa::SimpleTransposedOp<double>* isAtrans = dynamic_cast<const Playa::SimpleTransposedOp<double>* >(A.ptr().get());

	// Extract the underlying PlayaDenseSerialMatrix. Type checking is done
	// by rcp_dynamic_cast, so we need no error checking here.
	RCP<Playa::DenseSerialMatrix> APtr;
	if(isAtrans)
		APtr = rcp_dynamic_cast<DenseSerialMatrix>(isAtrans->op().ptr());	
	else
		APtr = rcp_dynamic_cast<DenseSerialMatrix>(A.ptr());
	
	// Get a pointer to the DenseSerialMatrix's data_
	const double* A0_ptr = APtr->dataPtr();

	// Scale the return value
	RCP<DenseSerialMatrix> MPtr;
	LinearOperator<double> rtnMtx;
	double* rtn;
	if(!isAtrans) // A was given
	{
		TEUCHOS_TEST_FOR_EXCEPTION(APtr->numCols() != d.dim(), std::runtime_error, "The " << APtr->numRows() << "x" << APtr->numCols() << " matrix A was given with the vector d with dimension " << d.dim() );
		MPtr = rcp(new DenseSerialMatrix(A.domain(), A.range()));
/*
 * This is necessary so that we can return our matrix as a LinearOperator.
 * Both LinearOperator and DenseSerialMatrix inherit from LinearOperatorBase
 */
		rtnMtx = rcp_dynamic_cast<LinearOperatorBase<double> >(MPtr);
		rtn = MPtr->dataPtr();
		for(int row = 0; row<APtr->numRows(); row++)
		{
			for(int col = 0; col<APtr->numCols(); col++)
			{
				rtn[row+APtr->numRows()*col] = A0_ptr[row+APtr->numRows()*col]*d[col];
			}
		}
	}
	else // A.transpose() was given
	{

		TEUCHOS_TEST_FOR_EXCEPTION(APtr->numRows() != d.dim(), std::runtime_error, "The " << APtr->numCols() << "x" << APtr->numRows() << " matrix A^T was given with the vector d with dimension " << d.dim() );
		MPtr = rcp(new DenseSerialMatrix(isAtrans->op().range(), isAtrans->op().domain()));
/*
 * This is necessary so that we can return our matrix as a LinearOperator.
 * Both LinearOperator and DenseSerialMatrix inherit from LinearOperatorBase
 */
		rtnMtx = rcp_dynamic_cast<LinearOperatorBase<double> >(MPtr);
		rtn = MPtr->dataPtr();
		for(int row = 0; row<APtr->numCols(); row++)
		{
			for(int col = 0; col<APtr->numRows(); col++)
			{
				rtn[row+APtr->numCols()*col] = A0_ptr[col+APtr->numRows()*row]*d[col];
			}
		}
	}	

 	return rtnMtx; 
}

/*
 * denseMatrixMatrixProduct performs C = A*B.
*/ 
LinearOperator<double> denseMatrixMatrixProduct(
  const LinearOperator<double>& A,
  const LinearOperator<double>& B)
{
// op(A) is m by kA, op(B) is kB by n and their product C is m by n
	char TRANSA = 'N';
	char TRANSB = 'N';
	int m = 1;
	int n = 1;
	int kA = 1;
	int kB = 0;
	RCP<Playa::DenseSerialMatrix> APtr;
	RCP<Playa::DenseSerialMatrix> BPtr;
	RCP<Playa::DenseSerialMatrix> CPtr;
	// means C = A*B
	double ALPHA = 1.0;
	double BETA = 0.0; 
	// Determine if either A or B are transposes
	const Playa::SimpleTransposedOp<double>* isAtrans = dynamic_cast<const Playa::SimpleTransposedOp<double>* >(A.ptr().get());
	const Playa::SimpleTransposedOp<double>* isBtrans = dynamic_cast<const Playa::SimpleTransposedOp<double>* >(B.ptr().get());
	if(isAtrans)
	{
		TRANSA = 'T';
		APtr = rcp_dynamic_cast<DenseSerialMatrix>(isAtrans->op().ptr());
		m = APtr->numCols();
		kA = APtr->numRows();
	}
	else
	{
		TRANSA = 'N';
		APtr = rcp_dynamic_cast<DenseSerialMatrix>(A.ptr());
		m = APtr->numRows();
		kA = APtr->numCols();
	}
	if(isBtrans)
	{
		TRANSB = 'T';
		BPtr = rcp_dynamic_cast<DenseSerialMatrix>(isBtrans->op().ptr());
		kB = BPtr->numCols();
		n = BPtr->numRows();
	}
	else
	{
		TRANSB = 'N';
		BPtr = rcp_dynamic_cast<DenseSerialMatrix>(B.ptr());
		kB = BPtr->numRows();
		n = BPtr->numCols();
	}
	

	TEUCHOS_TEST_FOR_EXCEPTION(kA != kB, std::runtime_error, "The dimensions of A and B are not compatible" );

	// Define the leading dimensions of A,B,C	
	int LDA = APtr->numRows();
	int LDB = BPtr->numRows();
	int LDC = m;

	//Define the dimensions of C
	if(isAtrans && isBtrans)
		CPtr = rcp(new DenseSerialMatrix(isBtrans->op().range(), isAtrans->op().domain()));
	else if(isAtrans)
		CPtr = rcp(new DenseSerialMatrix(B.domain(), isAtrans->op().domain()));
	else if(isBtrans)
		CPtr = rcp(new DenseSerialMatrix(isBtrans->op().range(), A.range()));
	else
		CPtr = rcp(new DenseSerialMatrix(B.domain(), A.range()));

	// dgemm performs C = alpha*A*B + beta*C
	dgemm_(&TRANSA, &TRANSB, &m, &n, &kA, &ALPHA, APtr->dataPtr(), &LDA, BPtr->dataPtr(), &LDB, &BETA, CPtr->dataPtr(), &LDC);


	LinearOperator<double> rtnMtx = rcp_dynamic_cast<LinearOperatorBase<double> >(CPtr);

 	return rtnMtx; 
}

}
