#include "PlayaSerialVectorType.hpp"
#include "PlayaDenseSerialMatrix.hpp"
#include "PlayaDenseLUSolver.hpp"
#include "Sundance.hpp"

//Needed to easily create A
#include <vector>
using std::vector;

//For using newton-armijo
#include "PlayaNewtonArmijoSolverImpl.hpp"

using std::cout;
using std::endl;

int main(int argc, char *argv[])
{

	int N = 3;
	VectorType<double> vecType = new SerialVectorType();
	VectorSpace<double> vecSpace = vecType.createEvenlyPartitionedSpace(MPIComm::self(), N);

	LinearOperator<double> A(rcp(new DenseSerialMatrix(vecSpace, vecSpace)));
	RCP<DenseSerialMatrix> APtr = DenseSerialMatrix::getConcretePtr(A);

	vector<vector<double> > AVals {{1,0,1},{2,1,2},{3,2,6}};
	APtr->fill(AVals);

	cout << "Here is A: " << endl << A << endl;

	Vector<double> b = vecSpace.createMember();
	b[0] = 3.0;
	b[1] = 7.0;
	b[2] = 8.0;

	Vector<double> x;

	DenseLUSolver solver;
	SolverState<double> state = solver.solve(A, b, x);
	TEUCHOS_TEST_FOR_EXCEPTION(state.finalState() != SolveConverged,
      	runtime_error, "solve failed");

    	Out::root() << "numerical solution = " << std::endl;
    	Out::os() << x << endl;

	string NLParamFile = "playa-newton-armijo.xml";
	ParameterXMLFileReader reader(NLParamFile);
	ParameterList solverParams = reader.getParameters();
	//Next line possible since DenseLUSolver and LinearSolver both inherit from LinearSolverBase
	LinearSolver<double> linearSolver(rcp(new DenseLUSolver())); 
	NewtonArmijoSolver<double> nonlinearSolver(solverParams, linearSolver);


	Vector<double> v = b;
	Vector<double> w = b.copy();
/*
	cout << "Here is b: " << endl << b << endl
             << "Here is v: " << endl << v << endl
             << "Here is w: " << endl << w << endl;
	b[1] = 100;
	cout << "Changed b[1] to 100" << endl;
	cout << "Here is b: " << endl << b << endl
             << "Here is v: " << endl << v << endl
             << "Here is w: " << endl << w << endl;
*/
	LinearOperator<double> B(rcp(new DenseSerialMatrix(vecSpace, vecSpace)));
	RCP<DenseSerialMatrix> Bptr = DenseSerialMatrix::getConcretePtr(B);
	double* BData = Bptr->dataPtr();
	for(int j=0; j<Bptr->numCols(); j++)
	{
		//Scale each column by j+1
		Vector<double> temp_b = (j+1.0)*b.copy();
		temp_b[1] = 2*j;	
		cout << "Value of temp_b for " << j << endl << temp_b << endl;
		temp_b*=(1.0/temp_b.norm2());
		cout << "Value of normalized temp_b for " << j << endl << temp_b;	
		cout << "2-norm value " << temp_b.norm2() << endl << endl;
		for(int i=0; i<Bptr->numRows(); i++)
			BData[i+Bptr->numRows()*j] = temp_b[i];
	}
	cout << "B is " << endl << B << endl;
/* Creating B with a vector did not tie the memory locations together
	cout << "Changing b[1] to 100" << endl;
	b[1] = 100;
	cout << "b is " << endl << b << endl;
	cout << "B is " << endl << B << endl;
*/
	cout << "Attempting to pull out each column of B" << endl;
	Vector<double> e_j = B.domain().createMember();
	Vector<double> b_j = B.range().createMember();
	for(int j=0; j<Bptr->numCols(); j++)
	{
		e_j.zero();
		e_j[j] = 1.0;
		b_j.zero();

		B.apply(e_j,b_j);
		cout << "Here is column " << j << " of B " << endl << b_j;
		cout << "2-norm value " << b_j.norm2() << endl << endl;
	}

	BasisArray basis = List(new Lagrange(2), new Lagrange(2));
	cout << "basis's size " << basis.size() << endl;



















return 0;
}
