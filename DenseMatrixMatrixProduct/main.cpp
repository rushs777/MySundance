#include "Teuchos_GlobalMPISession.hpp" //Handles initializing and finalizing MPI. Also defines rcp_dynamic_cast
#include "PlayaSerialVectorType.hpp" // For SerialVectorType
#include "PlayaDenseSerialMatrix.hpp" // For DenseSerialMatrix
#include "PlayaVectorType.hpp" // For VectorType, VectorSpace
#include "PlayaSimpleDiagonalOpDecl.hpp"// Needed for diagonalOperator
#include "PlayaSimpleComposedOpDecl.hpp"//Needed for matrix, matrix multiply via *
#include "PlayaLinearCombinationImpl.hpp"
//#include "PlayaOut.hpp" // What does this do?

#include "PlayaDenseMatrixMatrixProduct.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaSimpleTransposedOpDecl.hpp" // Needed if I want to take transposes.




#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp" // Needed to create an empty LinearOperator object
#include "PlayaVectorOpsImpl.hpp" // Needed for A*x
#include "PlayaVectorImpl.hpp" // Needed for A*x
#include "PlayaSimpleDiagonalOpImpl.hpp" // Needed for diagonalOperator
#include "PlayaSimpleComposedOpImpl.hpp" //Needed for matrix, matrix multiply via *
#include "PlayaSimpleTransposedOpImpl.hpp" //Needed if I want to take transposes
#endif



int main (int argc, char *argv[])
{      
	int stat = 0;
try
{

	Teuchos::GlobalMPISession session(&argc, &argv);


	int rows = 5;
	int cols = 3;

	Playa::VectorType<double> vecType = new Playa::SerialVectorType();
	Playa::VectorSpace<double> domain = vecType.createEvenlyPartitionedSpace(Playa::MPIComm::world(),cols);
	Playa::VectorSpace<double> range = vecType.createEvenlyPartitionedSpace(Playa::MPIComm::world(),rows);

	Teuchos::RCP<Playa::MatrixFactory<double> > mfA = vecType.createMatrixFactory(domain, range);
	Teuchos::RCP<Playa::MatrixFactory<double> > mfB = vecType.createMatrixFactory(range, domain);


	Playa::LinearOperator<double> A = mfA->createMatrix();
	Playa::LinearOperator<double> B = mfB->createMatrix();


      Teuchos::RCP<Playa::DenseSerialMatrix> PtrA	 
        = Teuchos::rcp_dynamic_cast<Playa::DenseSerialMatrix>(A.ptr());
      Teuchos::RCP<Playa::DenseSerialMatrix> PtrB	 
        = Teuchos::rcp_dynamic_cast<Playa::DenseSerialMatrix>(B.ptr());

	// Initialize test matrices
	// A will be 0 to rows*col-1, and B will be rows*cols-1
	for(int row = 0; row<PtrA->numRows(); row++)
		for(int col = 0; col<PtrA->numCols(); col++)
			PtrA->dataPtr()[row+PtrA->numRows()*col] = col+PtrA->numCols()*row;

	for(int row = 0; row<PtrB->numRows(); row++)
		for(int col = 0; col<PtrB->numCols(); col++)
			PtrB->dataPtr()[row+PtrB->numRows()*col] = col+PtrB->numCols()*row;
		
	std::cout << "Here's A: " << std::endl << A << std::endl
		  << "Here's B: " << std::endl << B << std::endl;



	// Transpose
	Playa::LinearOperator<double> At = A.transpose();
	Playa::LinearOperator<double> Bt = B.transpose();




	// Test LeftScale
	// testLeftScale will be 1 to rows
	Playa::Vector<double> testLeftScale = range.createMember();
	for(int i = 0; i<testLeftScale.dim(); i++)
	{
		testLeftScale[i] = i+1.0;
	}

	Playa::LinearOperator<double> resultLeftScale = denseLeftScale(testLeftScale,A);
	Playa::LinearOperator<double> resultLeftScalet = denseLeftScale(testLeftScale,Bt);

	// Test RightScale
	// testRightScale will be 1 to cols
	Playa::Vector<double> testRightScale = domain.createMember();
	for(int i = 0; i<testRightScale.dim(); i++)
	{
		testRightScale[i] = i+1.0;
	}
	Playa::LinearOperator<double> resultRightScale = denseRightScale(A,testRightScale);
	Playa::LinearOperator<double> resultRightScalet = denseRightScale(Bt,testRightScale);


	// Test denseMatrixMatrixProduct
	Playa::LinearOperator<double> resultDMMP = denseMatrixMatrixProduct(A,B);
	Playa::LinearOperator<double> resultDMtMP = denseMatrixMatrixProduct(At,A);
	Playa::LinearOperator<double> resultDMMtP = denseMatrixMatrixProduct(B,Bt);
	Playa::LinearOperator<double> resultDMtMtP = denseMatrixMatrixProduct(Bt,At);

	std::cout << "Here is AtA " << std::endl << resultDMtMP << std::endl;

	// Start of actual testing, using the ComposedOperator * to verify
	int nSamples = 10;
	bool allOK = true;
	double tol = 1.0e-12;
	for(int count=0; count<nSamples; count++)
	{
		std::cout << "Sample #" << count << " of " << nSamples << std::endl;
		Playa::Vector<double> x = domain.createMember();
		x.randomize();


		// Check denseLeftScale
		Playa::LinearOperator<double> TestLeftScale = Playa::diagonalOperator(testLeftScale);
		Playa::Vector<double> z = TestLeftScale*A*x - resultLeftScale*x;
		Playa::Vector<double> zt = TestLeftScale*Bt*x - resultLeftScalet*x;
		double errorz = z.norm2();
		double errorzt = zt.norm2();
		std::cout << "|| (D*A - denseLeftScale(d,A))*x || = " << errorz << std::endl;
		std::cout << "|| (D*Bt - denseLeftScale(d,Bt))*x || = " << errorzt << std::endl;

		// Check denseRightScale
		Playa::LinearOperator<double> TestRightScale = Playa::diagonalOperator(testRightScale);
		Playa::Vector<double> y = A*TestRightScale*x - resultRightScale*x;
		Playa::Vector<double> yt = Bt*TestRightScale*x - resultRightScalet*x;
		double errory = y.norm2();
		double erroryt = y.norm2();
		std::cout << "|| (A*D - denseRightScale(A,d))*x || = " << errory << std::endl;
		std::cout << "|| (Bt*D - denseRightScale(Bt,d))*x || = " << erroryt << std::endl;

		// Check denseMatrixMatrixProduct
		Playa::Vector<double> u = range.createMember();
		Playa::Vector<double> mm = A*B*u - resultDMMP*u;
		Playa::Vector<double> mtm = At*A*x - resultDMtMP*x;
		Playa::Vector<double> mmt = B*Bt*x - resultDMMtP*x;
		Playa::Vector<double> mtmt = Bt*At*u - resultDMtMtP*u;
		double errormm = mm.norm2();
		double errormtm = mtm.norm2();
		double errormmt = mmt.norm2();
		double errormtmt = mtmt.norm2();
		std::cout << "|| (A*B - denseMatrixMatrixProduct(A,B))*u || = " << errormm << std::endl;
		std::cout << "|| (At*A - denseMatrixMatrixProduct(At,A))*x || = " << errormtm << std::endl;
		std::cout << "|| (B*Bt - denseMatrixMatrixProduct(B,Bt))*x || = " << errormmt << std::endl;
		std::cout << "|| (Bt*At - denseMatrixMatrixProduct(Bt,At))*u || = " << errormtmt << std::endl;


		if (errorz > tol || errorzt > tol ||
		    errory > tol || erroryt > tol ||
		    errormm > tol || errormtm > tol || errormmt > tol || errormtmt > tol) 
			allOK = false;
	}

	      if (allOK)
	      {
	        Out::os() << "denseMatrixMatrixProduct test PASSED" << std::endl;
	      }
	      else
	      {
	        Out::os() << "denseMatrixMatrixProduct test FAILED" << std::endl;
	        stat = -1;
	      }


} //end of try

  catch(std::exception& e)
    {
      stat = -1;
      Out::os() << "Caught exception: " << e.what() << std::endl;
    }	



	return stat;
}
