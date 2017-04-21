#include "Teuchos_GlobalMPISession.hpp"
#include "QuadraticODERHSBase.hpp"

//For using newton-armijo
#include "PlayaDenseLUSolver.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "PlayaNewtonArmijoSolverImpl.hpp"
#include "PlayaNonlinearSolver.hpp"
#include "PlayaNonlinearSolverBuilder.hpp"
#include "MyNLO.hpp"

#include <vector>
using std::vector;

/*
  using namespace Teuchos;
  using namespace Playa;
  using namespace PlayaExprTemplates;
*/

class TestQuadODE : public QuadraticODERHSBase
{
public:
  TestQuadODE() : QuadraticODERHSBase(3, 2) {}

  Vector<double> evalForceTerm(const double& t) const
  {
    Vector<double> rtn = space().createMember();
    rtn[0] = 1 + t*(-7 + t*(5 + (-1 - t)*t));
    rtn[1] = 3 + t*(-3 + t*(-11 + t*(-5 + 2*t)));
    rtn[2] = -1 + t*(-4 + t*(6 + (-2 + t)*t));
    //rtn[0] = 1 + t*(-6 + t*(3 + (-1 - 3*t)*t));
    //rtn[1] = 4 + t*(4 + t*(-8 + t*(-2 + 2*t)));
    //rtn[2] = -2 + t*(-3 + t*(-1 + (-3 - 2*t)*t));
    return rtn;
  }

  void fillMatrixAndTensor(RCP<DenseSerialMatrix>& A,
			   Array<RCP<DenseSerialMatrix> >& T)
  {
    vector<vector<double> > AVals {{2,-1,0},{-1,2,-1},{0,-1,2}};
    vector<vector<vector<double> > > TVals = {{{-3,-4,4},{-3,1,0},{3,1,3}},
					      {{-2,3,4},{2,-2,3},{2,2,4}},
					      {{2,-3,4},{-1,-1,3},{-2,-1,-1}}};
    
    A->fill(AVals);
    for (int i=0; i<AVals.size(); i++)
      {
	T[i]->fill(TVals[i]);
      }
  }
};

/*****************************************************************************************
 *
 * The class MyNLO inherits from NonLinearOperatorBase. As such, it must provide
 * implementation for the pure virtual functions:
 *	Vector< double > getInitialGuess() const
 *	LinearOperator< double> computeJacobianAndFunction(Vector<double> &functionValue) const
 *
 * The purpose of MyNLO is work with the function
 *		F(z) = z - uPrev - (h/2.0)*(f(tPrev,uPrev) + f(tNext,z) )
 *
 * Which is the result of the Trapezoid Rule applied to a nonlinear ODE.
 * uNext = z + uPrev
 *
 *****************************************************************************************
class MyNLO : public NonlinearOperatorBase<double>
{
public:
  MyNLO(TestQuadODE f, double h) : NonlinearOperatorBase(f.space(), f.space()), f_(f) 
  {
    uPrev_ = f.space().createMember();
    setEvalPt(getInitialGuess());

    h_ = h;
    tPrev_ = 0.0;
    tNext_ = tPrev_+h_;
  }


  virtual RCP<NonlinearOperatorBase> getRcp() {return rcp(this);}

  Vector<double> getInitialGuess() const
  {
    return uPrev_.copy();
  }

  void set_tPrev(const double t)
  {
    tPrev_ = t;
    tNext_ = tPrev_ + h_;
  }

  void set_uPrev(const Vector<double> u)
  {
    uPrev_ = u.copy();
  }

protected:
  LinearOperator<double> computeJacobianAndFunction(Vector<double>& functionValue) const
  {
    LinearOperator<double> Jf;
    // Calculate f(t_n, u_n)
    Vector<double> fPrev = f_.eval1(tPrev_, uPrev_, Jf);
    // Calculate f(t_{n+1}, x^n)
    Vector<double> fNext = f_.eval1(tNext_, currentEvalPt(), Jf);
    // Calculate F(x^n)
    functionValue = currentEvalPt() - uPrev_ - (h_/2.0)*(fPrev + fNext);
    // Calculate JF(x^n) = I - (h/2)Jf(x^n)
    LinearOperator<double> JF(rcp(new DenseSerialMatrix(f_.space(), f_.space())));
    RCP<DenseSerialMatrix> JFptr = DenseSerialMatrix::getConcretePtr(JF);
    RCP<DenseSerialMatrix> Jfptr = DenseSerialMatrix::getConcretePtr(Jf);

    for(int i = 0; i<Jfptr->numRows(); i++)
      for(int j = 0; j<Jfptr->numCols(); j++)
	{
	  if(i==j)
	    JFptr->setElement(i,i, 1.0 - (h_/2.0)*Jfptr->getElement(i,i));
	  if(i!=j)
	    JFptr->setElement(i,j, (-1.0)*(h_/2.0)*Jfptr->getElement(i,j));
	}

    return JF;
  }

private:
  TestQuadODE f_;
  Vector<double> uPrev_;
  double h_;
  double tPrev_;
  double tNext_;
};
*/

int main(int argc, char *argv[]) 
{
  try
    {

      int verbosity = 0;	
      Sundance::setOption("verbosity", verbosity, "verbosity level");
      
      int nx = 32;
      Sundance::setOption("nx", nx, "Number of elements along each axis");

      int nSteps = 128;
      Sundance::setOption("nSteps", nSteps, "Number of time steps");

      double tFinal = 1.0;
      Sundance::setOption("tFinal", tFinal, "final time");

      Sundance::init(&argc, &argv);
      GlobalMPISession session(&argc, &argv);


      RCP<TestQuadODE> f = rcp(new TestQuadODE());
      f->initialize();
      
      string NLParamFile = "playa-newton-armijo.xml";
      ParameterXMLFileReader reader(NLParamFile);
      ParameterList solverParams = reader.getParameters();
      //Next line possible since DenseLUSolver and LinearSolver both inherit from LinearSolverBase
      LinearSolver<double> linearSolver(rcp(new DenseLUSolver())); 
      NewtonArmijoSolver<double> nonlinearSolver(solverParams, linearSolver);

      double deltat = tFinal/nSteps;
      MyNLO* prob = new MyNLO(f, deltat);
      NonlinearOperator<double> F = prob;
      //F.setVerb(verbosity);
      

      Array<Vector<double> > soln(nSteps+1);
      //Establish u(t_init)
      soln[0] = f->space().createMember();
      soln[0][0] = 1.0;
      soln[0][1] = 0.0;
      soln[0][2] = 0.0;
      prob->set_uPrev(soln[0]);
      
      for(int time = 1; time < soln.length(); time++)
	{
	  SUNDANCE_ROOT_MSG2(verbosity, "Nonlinear Solve for time step " + Teuchos::toString(time) + " of " + Teuchos::toString(nSteps));
	  prob->set_tPrev( (time-1.0)*deltat );
	  //Solve u(t_{time+1})
	  SolverState<double> state = nonlinearSolver.solve(F, soln[time]);
	  TEUCHOS_TEST_FOR_EXCEPTION(state.finalState() != SolveConverged,
				     runtime_error, "solve failed");
	  prob->set_uPrev(soln[time]);
	  if(verbosity>=2)
	    cout <<"Soln[" << time*deltat << "]" << endl << soln[time] << endl;
	}
      
      SUNDANCE_ROOT_MSG2(verbosity, "numerical solution finished");
      //Out::os() << soln[0] << std::endl;

      cout << "Compare numerical solution to exact solution " << endl;

      // Create uExact = {1, t^2,t}
      Vector<double> uExact = f->space().createMember();

      VectorType<double> serialVecType = new SerialVectorType();
      VectorSpace<double> serialVecSpace = serialVecType.createEvenlyPartitionedSpace(MPIComm::self(), nSteps+1);
      Vector<double> error = serialVecSpace.createMember();
     
      uExact[0] = 1.0;
      for(int time=0; time<error.dim(); time++)
	{
	  uExact[1] = (time*deltat)*(time*deltat);
	  uExact[2] = time*deltat;
	  error[time] = (uExact-soln[time]).norm2();
	  SUNDANCE_ROOT_MSG1(verbosity, "error[" + Teuchos::toString(time) + "] = " + Teuchos::toString(error[time]));
	}

      Out::root() << "||uExact - soln||_2  :\t " << error.norm2() << endl;
      Out::root() << "||uExact - soln||_inf:\t " << error.normInf() << endl;

      
    }
  catch(std::exception& e)
    {
      Out::root() << "Caught exception: " << e.what() << std::endl;
      return -1;
    }
}

