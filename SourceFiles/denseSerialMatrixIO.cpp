#include "denseSerialMatrixIO.hpp"

void writeDenseSerialMatrix(const RCP<DenseSerialMatrix>& A, string filename)
{
/* A needs to be a pointer to the underlying DenseSerialMatirx of a LinearOperator
   This is done via: LinearOperator<double> B(rcp(new DenseSerialMatrix(vecSpace, vecSpace)));	
	             RCP<DenseSerialMatrix> A = DenseSerialMatrix::getConcretePtr(B);
*/
ofstream os(filename);
	for(int i = 0; i<A->numRows(); i++)
		for(int j=0; j<A->numCols(); j++)
			os << A->getElement(i,j) << endl;
}


void readDenseSerialMatrix(RCP<DenseSerialMatrix>& A, string filename)
{
ifstream is(filename, std::ios::in);
	TEUCHOS_TEST_FOR_EXCEPTION(!is, std::runtime_error, "could not open file "
			     << filename << " for reading");

	double value;
	for(int i = 0; i<A->numRows(); i++)
		for(int j = 0; j<A->numCols(); j++)
		{
			is >> value;
			A->setElement(i,j,value);
		}
}


void readListExprMatrix(Expr &A, int N, string filename)
{
  ifstream is(filename, std::ios::in);
  TEUCHOS_TEST_FOR_EXCEPTION(!is, std::runtime_error, "could not open file " << filename << " for reading");

  // The information is written to file by denseSerialMatrixIO::writeDenseSerialMatrix which
  // is done row-wise
  double value;
  for(int i = 0; i < N; i++)
    {
      Expr row;
      for(int j = 0; j < N; j++)
	{
	  is >> value;
	  row.append(value);
	}
      A.append(row);
    }
}





Vector<double> meanVector(const LinearOperator<double> A)
{
  // Calculate the mean of the columns
  Vector<double> avgVec = A.range().createMember();
  Vector<double> ones = A.domain().createMember();
  ones.setToConstant(1.0);
  A.apply(ones, avgVec);
  avgVec *= (1.0/ A.domain().dim() );
  return avgVec;
}


Expr timeMeanFunctionGenerator(const LinearOperator<double> A, const DiscreteSpace ds)
{
  Expr avgFunc = new DiscreteFunction(ds, serialToEpetra(meanVector(A)));
  return avgFunc;
}





LinearOperator<double> generateFluctuationMatrix(const LinearOperator<double> A)
{
  Vector<double> avgVec = meanVector(A);

  Playa::VectorType<double> serialVecType = new Playa::SerialVectorType();
  RCP<MatrixFactory<double> > mf = serialVecType.createMatrixFactory(A.domain(), A.range());
  LinearOperator<double> rtnMtx = mf->createMatrix();
  RCP<DenseSerialMatrix> rtnMtxPtr = DenseSerialMatrix::getConcretePtr(rtnMtx);

  
  RCP<const DenseSerialMatrix> APtr = DenseSerialMatrix::getConcretePtr(A);
  double Aij;
  for(int j = 0; j < A.domain().dim(); j++)
    {
      for(int i = 0; i<A.range().dim(); i++)
	{
	  Aij = APtr->getElement(i,j);
	  rtnMtxPtr->setElement(i,j,Aij-avgVec[i]);
	}  
    }

  return rtnMtx;
}
