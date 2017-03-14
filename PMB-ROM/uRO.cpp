#include "Teuchos_GlobalMPISession.hpp"
#include "PlayaSerialVectorType.hpp"
#include "PlayaDenseSerialMatrix.hpp"
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaSerialEpetraAdapter.hpp" // Allows me to convert between Epetra and Serial Vectors
#include "PlayaOut.hpp"
#include <vector>


#include "VientoSnapshotIO.hpp" //For snapshotToMatrix
#include "denseSerialMatrixIO.hpp"
#include "QuadraticODERHSBase.hpp"

//For using newton-armijo
#include "PlayaDenseLUSolver.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "PlayaNewtonArmijoSolverImpl.hpp"
#include "PlayaNonlinearSolver.hpp"
#include "PlayaNonlinearSolverBuilder.hpp"


#include <sys/resource.h>
void memcheck()
{
  Tabs tab;
  struct rusage r_usage;
  getrusage(RUSAGE_SELF, &r_usage);
  Out::root() << tab << "Memory used (MB): " << r_usage.ru_maxrss / ((double) 1.0e6) << endl;
}


//Local files
#include "PlayaSVD.hpp"

#include "Sundance.hpp"

using std::vector;
using std::cout;
//using std::endl;
using namespace Teuchos;
using namespace Playa;
using namespace PlayaExprTemplates;

// Things necessary for mathematica generated fs to be understood
const double Pi = 4.0*atan(1.0);
Expr Cos(const Expr& x) {return cos(x);}
Expr Sin(const Expr& x) {return sin(x);}
Expr Power(const Expr& x, const double& p) {return pow(x,p);}
Expr Power(const double& x, const Expr& p) {return exp(p*log(x));}
const double E = exp(1.0);
Expr Sqrt(const Expr& x) {return sqrt(x);}
double Sqrt(const double& x) {return sqrt(x);}


class MMSQuadODE : public QuadraticODERHSBase
{
public:
  MMSQuadODE(Teuchos::Array<Expr> phi, Mesh mesh, bool MatrixAndTensorInFile = false, int verbosity = 1, int quadOrder = 6) 
    : QuadraticODERHSBase(phi.size(), verbosity),
      interior_(new MaximalCellFilter()),
      phi_(phi), mesh_(mesh),
      forceIP_(phi.size()),
      MatrixAndTensorInFile_(MatrixAndTensorInFile),
      quad_(new GaussianQuadrature(quadOrder))
  {
    t_ = new Sundance::Parameter(0.0);
    Expr x = new CoordExpr(0,"x");
    Expr y = new CoordExpr(1,"y");
    double nu = 1.0;
    q_ = List((-36*Cos(y)*Sin(t_)*Sin(x) + (((44 + 30*Cos(t_) + 8*Cos(4*t_) + 30*Cos(5*t_))*
					     Cos(x) + 9*((3 + 2*Cos(t_) + 2*Cos(5*t_))*Cos(3*x) + Cos(5*x)))*
					    Power(Sin(x),3) + 9*Cos(6*t_)*Cos(2*x)*Power(Sin(2*x),3))*Power(Sin(y),2) - 
	       3*(8*nu*Cos(2*t_)*(-1 + 2*Cos(2*x)) + 6*nu*Cos(3*t_)*(-1 + 5*Cos(4*x)) + 
		  8*Sin(2*t_)*Power(Sin(x),2) + 9*Sin(3*t_)*Power(Sin(2*x),2))*Sin(2*y))/36.,
	      (6*nu*(2*Cos(2*t_)*(-1 + 2*Cos(2*y))*Sin(2*x) + 
		     3*Cos(3*t_)*(-4 + 5*Cos(2*y))*Sin(4*x)) - 18*Cos(x)*Sin(t_)*Sin(y) + 
	       3*(4*Sin(2*t_)*Sin(2*x) + 9*Sin(3*t_)*Sin(4*x))*Power(Sin(y),2) + 
	       2*Cos(y)*(Cos(2*t_)*(4*Cos(2*t_) - 3*Cos(3*t_)*(-3 - 6*Cos(2*x) + Cos(4*x)))*
			 Power(Sin(x),2) + 9*Power(Cos(3*t_),2)*Power(Sin(2*x),2))*Power(Sin(y),3))/
	      18.);

    for(int r = 0; r < phi_.size(); r++)
      {
	forceIP_[r] = FunctionalEvaluator(mesh_, Integral(interior_, q_*phi_[r], quad_));
      }
  }

  Vector<double> evalForceTerm(const double& t) const
  {
    SUNDANCE_ROOT_MSG1(getVerbosity(), "start eval force");
    memcheck();
    t_.setParameterValue(t);

    Vector<double> rtn = space().createMember();
    SUNDANCE_ROOT_MSG1(getVerbosity(), "vec size: " << 8*space().dim());
    for(int r = 0; r < phi_.size(); r++)
      {
	rtn[r] = forceIP_[r].evaluate();
      }
    
    SUNDANCE_ROOT_MSG1(getVerbosity(), "end eval force");
    memcheck();
    //std::cout << "Here is the value of (vf, phi) " << std::endl << rtn << std::endl;
    return rtn;
  }

  void fillMatrixAndTensor(RCP<DenseSerialMatrix>& A, Array<RCP<DenseSerialMatrix> >& T)
  {
    Out::root() << "Starting fillMatrixAndTensor " << endl;
    // Access the matrix as a DenseSerialMatrix; note that A is a pointer to A_
    // This will create the matrix A = [(grad*phi_i, grad*phi_j)]
    string fileDir = "A_and_T";
    system( ("mkdir -p " + fileDir).c_str() );

    if(!MatrixAndTensorInFile_)
      {
	SUNDANCE_ROOT_MSG1(getVerbosity(), "Creating A");
	for(int i = 0; i<A->numRows(); i++)
	  for(int r = 0; r<A->numCols(); r++)
	    A->setElement(i,r,gradIP(phi_[i], phi_[r]));

	string A_filename = fileDir + "/A.txt";
	writeDenseSerialMatrix(A, A_filename);

	// Remember T is a pointer to T_
	for(int s = 0; s<phi_.size(); s++)
	  {
	    SUNDANCE_ROOT_MSG1(getVerbosity(), "Creating T[" + Teuchos::toString(s) + "] of " + Teuchos::toString(phi_.size()) );
	    for(int i = 0; i<A->numRows(); i++)
	      for(int r = 0; r<A->numCols(); r++)
		T[s]->setElement(i,r,tensorIP(phi_[i], phi_[s], phi_[r]));
	    string T_filename = fileDir + "/T[" + Teuchos::toString(s) + "].txt";		
	    writeDenseSerialMatrix(T[s], T_filename);
	  }
      }
    else
      {
	SUNDANCE_ROOT_MSG1(getVerbosity(), "Reading A from file");
	string A_filename = fileDir + "/A.txt";
	readDenseSerialMatrix(A, A_filename);

	for(int s = 0; s<phi_.size(); s++)
	  {
	    SUNDANCE_ROOT_MSG1(getVerbosity(), "Reading T[" + Teuchos::toString(s) + "] of " + Teuchos::toString(phi_.size()) + " from file");
	    string T_filename = fileDir + "/T[" + Teuchos::toString(s) + "].txt";	
	    readDenseSerialMatrix(T[s], T_filename);
	  }
      }
    SUNDANCE_ROOT_MSG1(getVerbosity(), "Done fillMatrixAndTensor");
  } // End of fillMatrixAndTensor
	

private:
  /********************************************************************************
   *
   * Teuchos::Array<Expr> phi - Array of basis function obtained from a POD;
   *                            Each Expr has an underlying DiscreteFunction
   * Mesh mesh - Considered the domain for the integral of the IP
   * QuadratureFamily quad - quadrature rule to use for the integral of the IP
   * mutable Expr t_ - time at this level is a Sundance::Paramter; mutable allows you to
   *                   affect t_ in functions that are labled as const
   * Expr q_ - Holds the volume force term
   * bool MatrixAndTensorInFile_ - True if A and T are on hand; false if they need to be created
   */
  CellFilter interior_;
  Teuchos::Array<Expr> phi_;
  Mesh mesh_;
  Teuchos::Array<FunctionalEvaluator> forceIP_;
  bool MatrixAndTensorInFile_;
  QuadratureFamily quad_;
  mutable Expr t_;
  Expr q_;

  /********************************************************************************
   * gradIP peforms the IP (grad*f, grad*g)
   ********************************************************************************/
  double gradIP(Expr f, Expr g)
  {
    // mesh.spatialDim() returns n for nD
    int dim = mesh_.spatialDim();

    // Define our differential operators; note Derivative(x=0)
    Expr grad = gradient(dim);

    FunctionalEvaluator IP = FunctionalEvaluator(mesh_, Integral(interior_, colonProduct(outerProduct(grad,f),outerProduct(grad,g)), quad_));
    return (IP.evaluate());
  }

  /********************************************************************************
   * tensorIP peforms the IP (f, (h*grad)*g)
   ********************************************************************************/
  double tensorIP(Expr f, Expr g, Expr h)
  {
    // mesh.spatialDim() returns n for nD
    // Define grad operator
    Expr grad = gradient(mesh_.spatialDim());

    //      The first one is what KL and I did on 1/27; then I realized the indices were off
    //	FunctionalEvaluator IP = FunctionalEvaluator(mesh, Integral(interior, h*((f*grad)*g),quad));
    FunctionalEvaluator IP = FunctionalEvaluator(mesh_, Integral(interior_, f*((h*grad)*g),quad_));
    return (IP.evaluate());
  }


}; // End of MMSQuadODE class



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
 *****************************************************************************************/
class MyNLO : public NonlinearOperatorBase<double>
{
public:
  //MyNLO(Teuchos::Array<Expr> phi, Mesh mesh) : NonlinearOperatorBase(),f_(phi,mesh) {}

  //MyNLO(Teuchos::Array<Expr> phi, Mesh mesh, const VectorSpace<double>& domain, const VectorSpace<double>& range) : NonlinearOperatorBase(domain, range), f_(phi,mesh) {}

  MyNLO(MMSQuadODE f, double h) : NonlinearOperatorBase(f.space(), f.space()), f_(f) 
  {
    uPrev_ = f.space().createMember();
    setEvalPt(getInitialGuess());

    h_ = h;
    tPrev_ = 0.0;
    tNext_ = h;
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
	    JFptr->setElement(i,j, (-h_/2.0)*Jfptr->getElement(i,j));
	}

    return JF;
  }

private:

  MMSQuadODE f_;
  Vector<double> uPrev_;
  Vector<double> uNext_;
  double h_;
  double tPrev_;
  double tNext_;
};


/*
 * 1. How to set up computeJacobainAndFunction so that the problem is only created once
 * 2. Clean up format/datatypes of ""
 * 3. How to calculate an actual initial guess
 * 4. How to handle updating the value for x^{n+1}
 *
 * Create a setuPrev and settPrev that are member functions of MyNLO but not of the public UI
 * NonlinearOperator; similar to lines 83 and 84 of PoissonBoltzmannTest.cpp
 * In ODERHS, make time a parameter for evalForceTerm. In MyNLO, make time a double
 */








int main(int argc, char *argv[]) 
{
  try
    {
      int verbosity = 1;	
      Sundance::setOption("verbosity", verbosity, "verbosity level");
      
      int nx = 64;
      Sundance::setOption("nx", nx, "Number of elements along each axis");

      int nSteps = 32;
      Sundance::setOption("nSteps", nSteps, "Number of time steps");

      double tFinal = 1.0;
      Sundance::setOption("tFinal", tFinal, "final time");

      bool AreMatrixAndTensorInFile = false;
      Sundance::setOption("MatrixAndTensorInFile", "MatrixAndTensorNotInFile", AreMatrixAndTensorInFile, "true if the matrix and tensor are available text files. false if they need to be created");

      Sundance::init(&argc, &argv);
      GlobalMPISession session(&argc, &argv);

      // Define our coordinates
      Expr x = new CoordExpr(0,"x");
      Expr y = new CoordExpr(1,"y");
      Expr t = new Sundance::Parameter(0.0);
      double deltat = tFinal/nSteps;

      // Define our solution
      Expr uExact = List((2*Cos(2*t)*Cos(y)*Power(Sin(x),2)*Sin(y))/3. + 
			 (Cos(3*t)*Cos(y)*Power(Sin(2*x),2)*Sin(y))/2.,
			 (-2*Cos(2*t)*Cos(x)*Sin(x)*Power(Sin(y),2))/3. - 
			 Cos(3*t)*Cos(2*x)*Sin(2*x)*Power(Sin(y),2));

      Expr pExact = Cos(x)*Cos(y)*Sin(t);

      // Define our mesh
      MeshType meshType = new Sundance::BasicSimplicialMeshType();
      double xmin = 0.0;
      double xmax = 2.0*Pi;
      double ymin = 0.0;
      double ymax = 2.0*Pi;
      MeshSource mesher = new Sundance::PartitionedRectangleMesher(xmin,xmax,nx,1,ymin,ymax,nx,1,meshType,0);
      Mesh mesh = mesher.getMesh();

      // Read the snapshots into a matrix
      //      string outDir = "/home/sirush/PhDResearch/ODETest/Results";
      string outDir = "../ODETest/Results";
      string fileDir = outDir + "/pcd-err-nx-" + Teuchos::toString(nx)
	+ "-nt-" + Teuchos::toString(nSteps);
      string tag = "st-v";
      string filename = fileDir + "/" + tag;
      Playa::LinearOperator<double> W = snapshotToMatrix(filename, nSteps, mesh);
      SUNDANCE_ROOT_MSG2(verbosity, "Size of W: " << W.range().dim() << " by " << W.domain().dim());

      // Perform the POD of the matrix W
      Playa::LinearOperator<double> U;
      Playa::LinearOperator<double> Phi;
      Playa::Vector<double> sigma;

      // Create a BasisFamily to express our solution
      Array<Sundance::BasisFamily> ubasis = List(new Sundance::Lagrange(2), new Sundance::Lagrange(2)); // 2nd order Piece-Wise Quad Lagrange in 2D
      // Define our vector type
      Playa::VectorType<double> epetraVecType = new Playa::EpetraVectorType();
      Sundance::DiscreteSpace ds(mesh, ubasis, epetraVecType);

      // W and mesh need to be defined
      SUNDANCE_ROOT_MSG1(verbosity, "Entering POD");
      POD(W,sigma,U,Phi,ds,verbosity);
      SUNDANCE_ROOT_MSG2(verbosity, "POD finished");

      Vector<double> lambda = sigma.copy();
      for(int count = 0; count < sigma.dim(); count++)
	lambda[count] = sqrt(lambda[count]);
      double lambdaTotal = lambda.norm1();
      double lambdaSum = 0.0;
      int R = 0;
      for(int i = 0; i<lambda.dim(); i++)
	{
	  lambdaSum += lambda[i];
	  if(lambdaSum/lambdaTotal >= .999)
	    {
	      // R is the number of lambdas to keep
	      R = i + 1;
	      SUNDANCE_ROOT_MSG2(verbosity, "Number of lambda[i] kept: " + Teuchos::toString(R));
	      break;
	    }
	}

      // Based off the value for R, create an appropriate VectorSpace<double>
      VectorType<double> R_vecType = new SerialVectorType();
      VectorSpace<double> R_vecSpace = R_vecType.createEvenlyPartitionedSpace(MPIComm::self(), R);




      // Looking at PlayaSVD.cpp, Phi is a DenseSerialMatrix
      Playa::Vector<double> ej = Phi.domain().createMember();
      Playa::Vector<double> phiCoeff = Phi.range().createMember(); // These are the coefficient vectors
      Array<Expr> phi(R); // These are the POD basis functions
      SUNDANCE_ROOT_MSG2(verbosity, "Size of phi: " + Teuchos::toString(phi.size()));

      //Needed for the integral
      CellFilter interior = new MaximalCellFilter();
      QuadratureFamily quad4 = new GaussianQuadrature(4);

      // Get the Expr phi_r(x)
      for(int r = 0; r<R; r++)
	{
	  ej.zero();
	  ej[r] = 1.0;
	  phiCoeff.zero();
	  Phi.apply(ej,phiCoeff);
	  phi[r] = new DiscreteFunction(ds, serialToEpetra(phiCoeff)); //DiscreteFunction requires Epetra vectors
	  TEUCHOS_TEST_FOR_EXCEPTION( fabs(L2Norm(mesh, interior, phi[r], quad4) - 1.0) >= 1.0e-6,
				     runtime_error, "||phi["+Teuchos::toString(r)+"]|| = " + Teuchos::toString(L2Norm(mesh, interior, phi[r], quad4)) + " != 1");
	}

      //Find the exact alphas
      Array<Vector<double> > alpha(nSteps+1);
      for(int count = 0; count<alpha.length(); count++)
	alpha[count] = R_vecSpace.createMember();

      // Needs to be of size nSteps+1;
      for(int tIndex=0; tIndex < alpha.length(); tIndex++)
	{
	  for(int r=0; r<R; r++)
	    {
	      // alpha_r(t_m) = ( uEx(t_m, x, y), phi[r] )
	      FunctionalEvaluator ExactEvaluator(mesh, Integral(interior, uExact*phi[r], quad4));
	      alpha[tIndex][r] = ExactEvaluator.evaluate();
	    }
	  t.setParameterValue(t.getParameterValue()+deltat);
	}

      
      SUNDANCE_ROOT_MSG1(verbosity, "Creating MMSQuadODE");
      // Create the nonlinear operator for solving our nonlinear ODE
      MMSQuadODE f(phi, mesh, AreMatrixAndTensorInFile, verbosity);
      f.initialize();

      
      SUNDANCE_ROOT_MSG1(verbosity, "Creating NLO");
      MyNLO* prob = new MyNLO(f, deltat);
      NonlinearOperator<double> F = prob;
      F.setVerb(verbosity);


      SUNDANCE_ROOT_MSG1(verbosity, "Creating solver");
      // create the Newton-Armijo solver
      string NLParamFile = "playa-newton-armijo.xml";
      ParameterXMLFileReader reader(NLParamFile);
      ParameterList params = reader.getParameters();
      const ParameterList& solverParams = params.sublist("NewtonArmijoSolver");
      //cout << "tau absolute " << solverParams.get<double>("Tau Absolute") << endl;
      cout << "the verbosity from playa...xm." << solverParams.get<int>("Verbosity") << endl;
	
      //Next line possible since DenseLUSolver and LinearSolver both inherit from LinearSolverBase
      LinearSolver<double> linearSolver(rcp(new DenseLUSolver())); 
      NewtonArmijoSolver<double> nonlinearSolver(solverParams, linearSolver);


      Array<Vector<double> > soln(nSteps+1);
      // Establish alpha(t_init)
      soln[0] = R_vecSpace.createMember();
      soln[0] = alpha[0].copy();
      prob->set_uPrev(soln[0]);
      
      SUNDANCE_ROOT_MSG1(verbosity, "Running...");      
      for(int time = 1; time < soln.length(); time++)
	{
	  SUNDANCE_ROOT_MSG1(verbosity, "Nonlinear Solve at time step " + Teuchos::toString(time) + " of " + Teuchos::toString(nSteps));
	  prob->set_tPrev( (time-1.0)*deltat );
	  SolverState<double> state = nonlinearSolver.solve(F, soln[time]);
	  TEUCHOS_TEST_FOR_EXCEPTION(state.finalState() != SolveConverged,
				     runtime_error, "solve failed");
	  prob->set_uPrev(soln[time]);
	}

      SUNDANCE_ROOT_MSG2(verbosity, "numerical solution finished");
      //Out::os() << soln[0] << std::endl;

      cout << "Compare numerical solution to exact solution " << endl;


      Vector<double> alphaError = Phi.domain().createMember();
      for(int i=0; i < alpha.length(); i++)
	{
	  cout << "Exact alpha(t=" << i << "): " << endl << alpha[i] << endl;	
	  cout << "Approximate alpha(t=" << i << "): " << endl << soln[i] << endl;
	  alphaError[i] = (alpha[i] - soln[i]).norm2();
	  cout << "Error for alpha(t=" << i << "): " << alphaError[i] << endl << endl;
	}
      cout << "Max error: " << alphaError.normInf() << endl;
      //cout << "lambda " << endl << lambda << endl;


      /*
      cout << "Staring to build reduced-order u from alphaExact " << endl;
      Array<Expr> uRO(nSteps+1);
      for(int n=0; n<alpha.length(); n++)
	{	
	  uRO[n] = List(0.0, 0.0);
	  for(int r=0; r<R; r++)
	    {
	      uRO[n] = uRO[n] + alpha[n][r]*phi[r]; 
	    }
	}

      cout << "Comparing uExact(t_n) to uRO(t_n)" << endl;
      Vector<double> l2norm = Phi.domain().createMember();
      cout << "size of l2norm: " << l2norm.dim() << endl;
      for(int time = 0; time < uRO.length(); time++)
	{
	  t.setParameterValue(time*deltat);
	  l2norm[time] = L2Norm(mesh, interior, uExact - uRO[time], quad4);
	  cout << "Error at time " << time*deltat << "= " << l2norm[time] << endl;
	}
      */

      
      /*
      // Visualize the results
      string vtkDir = "Results/Visuals/";
      string vtkfilename = "nx"+Teuchos::toString(nx)+"nt"+Teuchos::toString(nSteps);
      system( ("mkdir -p " + vtkDir).c_str() ); 
      FieldWriter writer = new VTKWriter(vtkDir+vtkfilename);
      writer.addMesh(mesh);

      L2Projector projector(ds, uExact);
      //Write uExact for all the time steps
      for(int time=0; time<1; time++)
	{
	  t.setParameterValue(time*deltat);
	  L2Projector projectorRO(ds, uRO[time]);
	  writer.addField("uExact[0]", new ExprFieldWrapper(projector.project()[0]) );
	  writer.addField("uExact[1]", new ExprFieldWrapper(projector.project()[1]) );
	  writer.addField("uRO[0]", new ExprFieldWrapper(projectorRO.project()[0]) );
	  writer.addField("uRO[1]", new ExprFieldWrapper(projectorRO.project()[1]) );
	  //writer.addField("error", new ExprFieldWrapper(l2norm[time]) );
	  writer.write();
	}
      */

      //      cout << "The 2-norm for the velocity error over all " << nSteps+1 << " timesteps: " << l2norm.norm2() << endl;
	
    }
  catch(std::exception& e)
    {
      Out::root() << "Caught exception: " << e.what() << std::endl;
      return -1;
    }
}

