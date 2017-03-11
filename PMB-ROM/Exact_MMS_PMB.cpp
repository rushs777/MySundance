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


//Local files
#include "PlayaSVD.hpp"
#include "MathematicaConverter.hpp" // converts Mathematica ouput to c++

#include "Sundance.hpp"

using std::vector;
using std::cout;
//using std::endl;
using namespace Teuchos;
using namespace Playa;
using namespace PlayaExprTemplates;

class MMSQuadODE : public QuadraticODERHSBase
{
public:
  MMSQuadODE(Teuchos::Array<Expr> phi, Mesh mesh, bool MatrixAndTensorInFile = false, int verbosity = 1, int quadOrder = 6) 
    : QuadraticODERHSBase(phi.size(), verbosity), 
      phi_(phi), mesh_(mesh), 
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
  }

  Vector<double> evalForceTerm(const double& t) const
  {
    t_.setParameterValue(t);
    CellFilter interior = new MaximalCellFilter();

    Vector<double> rtn = space().createMember();
    for(int r = 0; r < phi_.size(); r++)
      {
	FunctionalEvaluator IP = FunctionalEvaluator(mesh_, Integral(interior, q_*phi_[r], quad_));		
	rtn[r] = IP.evaluate();
      }

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
  Teuchos::Array<Expr> phi_; 
  Mesh mesh_;
  bool MatrixAndTensorInFile_;
  QuadratureFamily quad_;
  mutable Expr t_;
  Expr q_;

  /********************************************************************************
   * gradIP peforms the IP (grad*f, grad*g)
   ********************************************************************************/
  double gradIP(Expr f, Expr g)
  {

    CellFilter interior = new MaximalCellFilter();
    // mesh.spatialDim() returns n for nD
    int dim = mesh_.spatialDim();

    // Define our differential operators; note Derivative(x=0)
    Expr grad = gradient(dim);

    FunctionalEvaluator IP = FunctionalEvaluator(mesh_, Integral(interior, colonProduct(outerProduct(grad,f),outerProduct(grad,g)), quad_));
    return (IP.evaluate());
  }

  /********************************************************************************
   * tensorIP peforms the IP (f, (h*grad)*g)
   ********************************************************************************/
  double tensorIP(Expr f, Expr g, Expr h)
  {
    CellFilter interior = new MaximalCellFilter();
    // mesh.spatialDim() returns n for nD
    // Define grad operator
    Expr grad = gradient(mesh_.spatialDim());

    //      The first one is what KL and I did on 1/27; then I realized the indices were off
    //	FunctionalEvaluator IP = FunctionalEvaluator(mesh, Integral(interior, h*((f*grad)*g),quad));
    FunctionalEvaluator IP = FunctionalEvaluator(mesh_, Integral(interior, f*((h*grad)*g),quad_));
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



// Needed for the peg value; this isolates (0,0)
CELL_PREDICATE(CornerPointTest, {return fabs(x[1]) < 1.0e-10 && fabs(x[0])<1.0e-10;})




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
      string outDir = "/home/sirush/PhDResearch/ODETest/Results";
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
	  TEUCHOS_TEST_FOR_EXCEPTION(L2Norm(mesh, interior, phi[r], quad4) == 1.0,
				     runtime_error, "||phi["+Teuchos::toString(r)+"]|| != 1");
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
      for(int time = 0; time < uRO.length(); time++)
	{
	  t.setParameterValue(time*deltat);
	  l2norm[time] = L2Norm(mesh, interior, uExact - uRO[time], quad4);
	  cout << "Error at time " << time*deltat << "= " << l2norm[time] << endl;
	}

      double nu = 1.0;
      Expr q = List((-36*Cos(y)*Sin(t)*Sin(x) + (((44 + 30*Cos(t) + 8*Cos(4*t) + 30*Cos(5*t))*
						   Cos(x) + 9*((3 + 2*Cos(t) + 2*Cos(5*t))*Cos(3*x) + Cos(5*x)))*Power(Sin(x),3) + 9*Cos(6*t)*Cos(2*x)*Power(Sin(2*x),3))*Power(Sin(y),2) - 3*(8*nu*Cos(2*t)*(-1 + 2*Cos(2*x)) + 6*nu*Cos(3*t)*(-1 + 5*Cos(4*x)) + 8*Sin(2*t)*Power(Sin(x),2) + 9*Sin(3*t)*Power(Sin(2*x),2))*Sin(2*y))/36., (6*nu*(2*Cos(2*t)*(-1 + 2*Cos(2*y))*Sin(2*x) + 3*Cos(3*t)*(-4 + 5*Cos(2*y))*Sin(4*x)) - 18*Cos(x)*Sin(t)*Sin(y) + 3*(4*Sin(2*t)*Sin(2*x) + 9*Sin(3*t)*Sin(4*x))*Power(Sin(y),2) + 2*Cos(y)*(Cos(2*t)*(4*Cos(2*t) - 3*Cos(3*t)*(-3 - 6*Cos(2*x) + Cos(4*x)))*
Power(Sin(x),2) + 9*Power(Cos(3*t),2)*Power(Sin(2*x),2))*Power(Sin(y),3))/18.);
      
      // BasisFamily for p
      //BasisArray pbasis;
      //pbasis.push_back(new Sundance::Lagrange(1));
      BasisFamily pbasis = new Sundance::Lagrange(1);
      
      // Test Function for p
      //Sundance::Expr v = new Sundance::TestFunction(pbasis, "v");

      // Unknown Function for p
      //Sundance::Expr p = new Sundance::UnknownFunction(pbasis, "p");

      // mesh.spatialDim() returns n for nD
      int dim = mesh.spatialDim();

      // Define our differential operators; note Derivative(x=0)
      Expr grad = gradient(dim);

      /*******************Begin trying to do PMB****************************************/

      string pressureTag = "st-p";
      string pressureFilename = fileDir + "/" + pressureTag;
      Playa::LinearOperator<double> Q = snapshotToMatrix(pressureFilename, nSteps, mesh);

      // Perform the POD of the matrix Q
      Playa::LinearOperator<double> V;
      Playa::LinearOperator<double> Psi;
      Playa::Vector<double> pressureSigma;
      DiscreteSpace p_ds(mesh, pbasis, epetraVecType);

      // W and mesh need to be defined
      SUNDANCE_ROOT_MSG1(verbosity, "Entering POD for pressure");
      POD(Q,pressureSigma,V,Psi,p_ds,verbosity);
      SUNDANCE_ROOT_MSG2(verbosity, "POD finished for pressure");

      Vector<double> pressureLambda = pressureSigma.copy();
      for(int count = 0; count < pressureSigma.dim(); count++)
	pressureLambda[count] = sqrt(pressureLambda[count]);
      lambdaTotal = pressureLambda.norm1();
      lambdaSum = 0.0;
      int pressure_R = 0;
      for(int i = 0; i<pressureLambda.dim(); i++)
	{
	  lambdaSum += pressureLambda[i];
	  if(lambdaSum/lambdaTotal >= .999)
	    {
	      // R is the number of lambdas to keep
	      pressure_R = i + 1;
	      SUNDANCE_ROOT_MSG2(verbosity, "Number of pressureLambda[i] kept: " + Teuchos::toString(pressure_R));
	      break;
	    }
	}

      // Based off the value for R, create an appropriate VectorSpace<double>
      VectorType<double> pressure_R_vecType = new SerialVectorType();
      VectorSpace<double> pressure_R_vecSpace = pressure_R_vecType.createEvenlyPartitionedSpace(MPIComm::self(), pressure_R);

      // Looking at PlayaSVD.cpp, Psi is a DenseSerialMatrix
      Playa::Vector<double> ek = Psi.domain().createMember();
      Playa::Vector<double> psiCoeff = Psi.range().createMember(); // These are the coefficient vectors
      Array<Expr> psi(pressure_R); // These are the pressure POD basis functions
      SUNDANCE_ROOT_MSG2(verbosity, "Size of psi: " + Teuchos::toString(psi.size()));

      // Get the Expr psi_r(x)
      for(int r = 0; r<pressure_R; r++)
	{
	  ek.zero();
	  ek[r] = 1.0;
	  psiCoeff.zero();
	  Psi.apply(ek,psiCoeff);
	  psi[r] = new DiscreteFunction(p_ds, serialToEpetra(psiCoeff)); //DiscreteFunction requires Epetra vectors
	  TEUCHOS_TEST_FOR_EXCEPTION(L2Norm(mesh, interior, psi[r], quad4) == 1.0,
				     runtime_error, "||psi["+Teuchos::toString(r)+"]|| != 1");
	}

      //Find the exact betas
      Array<Vector<double> > beta(nSteps+1);
      for(int count = 0; count<beta.length(); count++)
	beta[count] = pressure_R_vecSpace.createMember();

      t.setParameterValue(0.0);
      for(int tIndex=0; tIndex < beta.length(); tIndex++)
	{
	  for(int r=0; r<pressure_R; r++)
	    {
	      // beta_r(t_m) = ( pEx(t_m, x, y), psi[r] )
	      FunctionalEvaluator ExactEvaluator(mesh, Integral(interior, pExact*psi[r], quad4));
	      beta[tIndex][r] = ExactEvaluator.evaluate();
	    }
	  t.setParameterValue(t.getParameterValue()+deltat);
	}
      

      cout << "Staring to build reduced-order p from betaExact " << endl;
      Array<Expr> pRO(nSteps+1);
      for(int m=0; m<beta.length(); m++)
	{	
	  pRO[m] = 0.0;
	  for(int r=0; r<pressure_R; r++)
	    {
	      pRO[m] = pRO[m] + beta[m][r]*psi[r]; 
	    }
	}

      cout << "Comparing pExact(t_n) to pRO(t_n)" << endl;
      Vector<double> pressure_l2norm = Psi.domain().createMember();
      for(int time = 0; time < pRO.length(); time++)
	{
	  t.setParameterValue(time*deltat);
	  pressure_l2norm[time] = L2Norm(mesh, interior, pExact - pRO[time], quad4);
	  cout << "Error at time " << time*deltat << "= " << pressure_l2norm[time] << endl;
	}
 
      /*
           Vector<double> l2norm_p = Phi.domain().createMember();
      for(int time = 0; time < pRO.length(); time++)
	{
	  t.setParameterValue(time*deltat);
	  l2norm_p[time] = L2Norm(mesh, interior, pExact - pRO[time], quad4);
	  cout << "Error at time " << time*deltat << "= " << l2norm_p[time] << endl;
	}

      // Visualize the results
      string vtkDir = "Results/Visuals/";
      system( ("mkdir -p " + vtkDir).c_str() ); 
      string vtkfilename = "nx"+Teuchos::toString(nx)+"nt"+Teuchos::toString(nSteps);
      FieldWriter writer = new VTKWriter(vtkDir+vtkfilename);
      writer.addMesh(mesh);

      L2Projector uEx_projector(ds, uExact);

      L2Projector pEx_projector(p_ds, pExact);
      //Write uExact for all the time steps
      for(int time=0; time<1; time++)
	{
	  t.setParameterValue(time*deltat);
	  L2Projector uRO_projector(ds, uRO[time]);
	  L2Projector pRO_projector(p_ds, pRO[time]);
	  L2Projector uErrorProjector(ds, uExact - uRO[time]);
	  L2Projector pErrorProjector(p_ds, pExact - pRO[time]);
	  writer.addField("uExact[0]", new ExprFieldWrapper(uEx_projector.project()[0]) );
	  writer.addField("uExact[1]", new ExprFieldWrapper(uEx_projector.project()[1]) );
	  writer.addField("pExact", new ExprFieldWrapper(pEx_projector.project()) );
	  writer.addField("uRO[0]", new ExprFieldWrapper(uRO_projector.project()[0]) );
	  writer.addField("uRO[1]", new ExprFieldWrapper(uRO_projector.project()[1]) );
	  writer.addField("pRO", new ExprFieldWrapper(pRO_projector.project()) );
	  writer.addField("uError[0]", new ExprFieldWrapper(uErrorProjector.project()[0]) );
	  writer.addField("uError[1]", new ExprFieldWrapper(uErrorProjector.project()[1]) );
	  writer.addField("pError", new ExprFieldWrapper(pErrorProjector.project()) );
	  writer.write();
	}
      */

      cout << "The 2-norm for the velocity error over all timesteps: " << l2norm.norm2() << endl;
      cout << "The 2-norm for the pressure error over all timestpes: " << pressure_l2norm.norm2() << endl;
	
    }
  catch(std::exception& e)
    {
      Out::root() << "Caught exception: " << e.what() << std::endl;
      return -1;
    }
}

