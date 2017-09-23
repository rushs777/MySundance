#include "Sundance.hpp"
#include "PlayaDenseLUSolver.hpp"

// Local Files
#include "MyNLO.hpp"
#include "QuadraticODERHSBase.hpp"

// Standard Library Functions
using std::vector;
using std::cout;
using std::endl;

// Namespaces
using namespace Teuchos;
using namespace Playa;

/**
 * The class ODETest_1D is meant to implement the RHS of the ODE x'(t) = x(t)
 * A = 1, T = 0, b = 0 (force term) 
 */
class ODETest_1D : public QuadraticODERHSBase
{
public:
  ODETest_1D() : QuadraticODERHSBase(1,2) {} // (size of A, verbosity)

  // evalForceTerm is a pure void function in QuadraticODERHSBase
  Vector<double> evalForceTerm(const double& t) const
  {
    Vector<double> rtn = space().createMember();
    rtn[0] = 0.0;

    return rtn;
  }

  // fillMatrixAndTensor is a pure void function in QuadraticODERHSBase
  void fillMatrixAndTensor(RCP<DenseSerialMatrix>& A,
			   Array<RCP<DenseSerialMatrix> >& T)
  {
    vector<vector<double> > AVals {{1}};
    vector<vector<vector<double> > > TVals {{{0}}};

    A->fill(AVals);
    for(int i=0; i < AVals.size(); i++)
      {
	T[i]->fill(TVals[i]);
      }
    
  }
  
};

/**
 * The class ODETest2_1D is meant to implement the RHS of the 
 * ODE x'(t) = x(t) - (x(t))^2
 * A = 1, T = -1, b = 0 (force term) 
 */
class ODETest2_1D : public QuadraticODERHSBase
{
public:
  ODETest2_1D() : QuadraticODERHSBase(1,2) {} // (size of A, verbosity)

  // evalForceTerm is a pure void function in QuadraticODERHSBase
  Vector<double> evalForceTerm(const double& t) const
  {
    Vector<double> rtn = space().createMember();
    rtn[0] = 0.0;

    return rtn;
  }

  // fillMatrixAndTensor is a pure void function in QuadraticODERHSBase
  void fillMatrixAndTensor(RCP<DenseSerialMatrix>& A,
			   Array<RCP<DenseSerialMatrix> >& T)
  {
    vector<vector<double> > AVals {{1}};
    vector<vector<vector<double> > > TVals {{{-1}}};

    A->fill(AVals);
    for(int i=0; i < AVals.size(); i++)
      {
	T[i]->fill(TVals[i]);
      }
    
  }
  
};

int main(int argc, char *argv[])
{
  try
    {
      // Command line options
      // Can't give x0 Parameter values, so give it a double for now
      double alpha = 1.0;
      Sundance::setOption("alpha",alpha,"Initial value for x(0)");

      int verbosity = 1;
      Sundance::setOption("verbosity",verbosity,"Level of verbosity for displaying output");

      double tInit = 0.0;
      Sundance::setOption("tInit", tInit, "Initial time for the ODE");

      double T = 2.0;
      Sundance::setOption("T",T,"Final time for the ODE");

      int nSteps = 10;
      Sundance::setOption("nSteps",nSteps,"Number of time steps to go from tInit to T");
      
      // Initialize Sundance
      Sundance::init(&argc, &argv);

      RCP<ODETest_1D> f = rcp(new ODETest_1D());
      f->initialize();

      // Create the inital condition x(tInit) = alpha
      Vector<double> x0 = f->space().createMember();
      x0[0] = alpha;

      // Determine the size of deltaT
      double deltaT = (T - tInit)/nSteps;

      // Displays as position tab value
      //cout << "Here is x0: " << x0 << endl;

      SUNDANCE_ROOT_MSG1(verbosity, "Creating the nonlinear operator");
      MyNLO* prob = new MyNLO(f, deltaT);
      NonlinearOperator<double> F = prob; // This is F(w) = 0

      SUNDANCE_ROOT_MSG1(verbosity, "Creating nonlinear solver");
      ParameterXMLFileReader reader("playa-newton-armijo.xml");
      ParameterList params = reader.getParameters();
      const ParameterList& nlSolverParams = params.sublist("NewtonArmijoSolver");
      //Next line possible since DenseLUSolver and LinearSolver both inherit from LinearSolverBase
      LinearSolver<double> linearSolver(rcp(new DenseLUSolver()));
      NewtonArmijoSolver<double> nonlinearSolver(nlSolverParams, linearSolver);

      SUNDANCE_ROOT_MSG1(verbosity,"Running nonlinear solves...");
      Array<Vector<double> > xSoln(nSteps+1);
      xSoln[0] = x0; // Only need to give the initial; solve expects empty Vec
      
      prob->set_uPrev(x0);
      for(int time = 1; time <= nSteps; time++)
	{
	  // Find x(t_{time})
	  SUNDANCE_ROOT_MSG1(verbosity, "Nonlinear Solver at time step " + Teuchos::toString(time) + " of " + Teuchos::toString(nSteps));
	  prob->set_tPrev( tInit + (time-1.0)*deltaT );
	  SolverState<double> state = nonlinearSolver.solve(F,xSoln[time]);

	  TEUCHOS_TEST_FOR_EXCEPTION(state.finalState() != SolveConverged,
				     runtime_error, "solve failed");

	  // Store x_{time} as the value for uPrev
	  prob->set_uPrev(xSoln[time]);	  
	}

      // xExact
      Expr t = new Sundance::Parameter(tInit);
      Expr xExact = alpha*exp(t);

      // Check xSoln against xExact
      VectorType<double> time_vecType = new SerialVectorType();
      VectorSpace<double> time_vecSpace = time_vecType.createEvenlyPartitionedSpace(MPIComm::self(), nSteps+1.0);

      Vector<double> xError = time_vecSpace.createMember();
      Vector<double> tempExact = x0.copy();
      for(int i=0; i <= nSteps; i++)
	{
	  /* fix this
	  t.setParameterValue(tInit + i*deltaT);
	  xError[i] = (xExact - xSoln[i]).norm2();
	  */
	  tempExact[0] = alpha*exp(tInit + i*deltaT);
	  xError[i] = (tempExact - xSoln[i]).norm2();
	  if(verbosity>=2)
	    cout << "Error for x(" << tInit + i*deltaT << "): " << xError[i] << endl;
	  cout << "xExact(" << tInit + i*deltaT << "): " << tempExact[0]
	       << "\t xSoln("  << tInit + i*deltaT << "): " << xSoln[i][0] << endl;
	}

      cout << "||vector of 2norms of the error in x(t_n), n:0:nSteps=" << nSteps << "||_2: " << xError.norm2() << endl;










      /*************************************************************
       *
       * Important note about this problem. An alpha value of 1.0
       * causes -alpha + 1 in the denominator to cancel out, and 
       * at that point we have have ( alpha*Exp(t) ) / (alpha*Exp(t))
       * which is just 1. So pick a different value for alpha to check
       * code
       *
       ************************************************************/

      RCP<ODETest2_1D> f2 = rcp(new ODETest2_1D());
      f2->initialize();

      // Use the same initial conditions x(tInit) = alpha
      // Use the same deltaT
      
      SUNDANCE_ROOT_MSG1(verbosity, "Creating the nonlinear operator");
      MyNLO* prob2 = new MyNLO(f2, deltaT);
      NonlinearOperator<double> F2 = prob2; // This is F(w) = 0

      SUNDANCE_ROOT_MSG1(verbosity, "Creating nonlinear solver");
      // Use the same nonlinear solver parameters
      SUNDANCE_ROOT_MSG1(verbosity,"Running nonlinear solves...");
      Array<Vector<double> > xSoln2(nSteps+1);
      xSoln2[0] = x0; // Only need to give the initial; solve expects empty Vec
      
      prob2->set_uPrev(x0);
      for(int time = 1; time <= nSteps; time++)
	{
	  // Find x(t_{time})
	  SUNDANCE_ROOT_MSG1(verbosity, "Nonlinear Solver at time step " + Teuchos::toString(time) + " of " + Teuchos::toString(nSteps));
	  prob2->set_tPrev( tInit + (time-1.0)*deltaT );
	  SolverState<double> state2 = nonlinearSolver.solve(F2,xSoln2[time]);

	  TEUCHOS_TEST_FOR_EXCEPTION(state2.finalState() != SolveConverged,
				     runtime_error, "solve failed");

	  // Store x_{time} as the value for uPrev
	  prob2->set_uPrev(xSoln2[time]);	  
	}

      // xExact2
      t.setParameterValue(tInit);
      Expr xExact2 = (alpha * exp(t) ) / (-alpha + alpha*exp(t) + 1.0);

      // Check xSoln2 against xExact2
      Vector<double> xError2 = time_vecSpace.createMember();
      Vector<double> tempExact2 = x0.copy();
      for(int i=0; i <= nSteps; i++)
	{
	  /* fix this
	  t.setParameterValue(tInit + i*deltaT);
	  xError[i] = (xExact - xSoln[i]).norm2();
	  */
	  tempExact2[0] = ( alpha*exp(tInit + i*deltaT) ) / ( -alpha + alpha*exp(tInit + i*deltaT) + 1.0);
	  xError2[i] = (tempExact2 - xSoln2[i]).norm2();
	  if(verbosity>=2)
	    cout << "Error for x(" << tInit + i*deltaT << "): " << xError2[i] << endl;
	  cout << "xExact(" << tInit + i*deltaT << "):\t " << tempExact2[0]
	       << "\t xSoln("  << tInit + i*deltaT << "):\t " << xSoln2[i][0] << endl;
	}

      cout << "||vector of 2norms of the error in x(t_n), n:0:nSteps=" << nSteps << "||_2: " << xError2.norm2() << endl;
      
      Sundance::finalize();
      
    }
  catch(std::exception& e)
    {
      Out::root() << "Caught exception: " << e.what() << endl;
      return -1;
    }
}
