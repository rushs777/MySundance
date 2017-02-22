#include "denseSerialMatrixIO.hpp"

using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
	Playa::VectorType<double> vecType = new Playa::SerialVectorType();
	VectorSpace<double> vecSpace = vecType.createEvenlyPartitionedSpace(MPIComm::self(),3);
	LinearOperator<double> A( rcp(new DenseSerialMatrix(vecSpace, vecSpace) ) );	
	RCP<DenseSerialMatrix> APtr = DenseSerialMatrix::getConcretePtr(A);	
	vector<vector<double> > AVals {{2,-1,0},{-1,2,-1},{0,-1,2}};
	APtr->fill(AVals);

	cout << "Here is A: " << endl << A << endl;

	string filename = "matrix_test.txt";
	writeDenseSerialMatrix(APtr, filename);


	LinearOperator<double> B( rcp(new DenseSerialMatrix(vecSpace, vecSpace) ) );
	RCP<DenseSerialMatrix> BPtr = DenseSerialMatrix::getConcretePtr(B);	
	readDenseSerialMatrix(BPtr, filename);
	cout << "Here is B " << endl << B << endl;
	

	return 0;
}
