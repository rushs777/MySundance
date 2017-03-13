#include "denseSerialMatrixIO.hpp"
/** */
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

/** */
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


