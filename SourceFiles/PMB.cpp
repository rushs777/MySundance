#include "PMB.hpp"

PMB::PMB(const string filename, const Array<Expr>& uRO, Expr forceTerm, Expr t, const DiscreteSpace& ds, double deltat, double tolerance, int verbosity) : filename_(filename), uRO_(uRO), forceTerm_(forceTerm), t_(t), ds_(ds), deltat_(deltat), M_(uRO.length()), tol_(tolerance), verbosity_(verbosity)
{
  t_.setParameterValue(0.0);
}

void PMB::initialize()
{
  int Mtemp = M_ - 1;
  Qprime_ = snapshotToMatrix(filename_, Mtemp, ds_.mesh());

  // Calculate pbar(x)
  Vector<double> pbarVec = Qprime_.range().createMember();
  Vector<double> ones = Qprime_.domain().createMember();
  ones.setToConstant(1.0);
  Qprime_.apply(ones, pbarVec);
  pbarVec *= (1.0/M_);
  pbar_ = new DiscreteFunction(ds_, serialToEpetra(pbarVec));

  // Make Q be Q'
  RCP<DenseSerialMatrix> QPtr = DenseSerialMatrix::getConcretePtr(Qprime_);
  double Qij;
  for(int j = 0; j < Qprime_.domain().dim(); j++)
    {
      for(int i = 0; i<Qprime_.range().dim(); i++)
	{
	  Qij = QPtr->getElement(i,j);
	  QPtr->setElement(i,j,Qij-pbarVec[i]);
	}  
    }
      
      // Perform the POD of the matrix Q'
      Playa::LinearOperator<double> V;
      Playa::LinearOperator<double> Psi;
      Playa::Vector<double> sigma;


      // Q' and mesh need to be defined
      SUNDANCE_ROOT_MSG1(verbosity_, "Entering POD for pressure");
      POD(Qprime_,sigma,V,Psi,ds_,verbosity_);
      SUNDANCE_ROOT_MSG2(verbosity_, "POD finished for pressure");

      Vector<double> lambda = sigma.copy();
      for(int count = 0; count < sigma.dim(); count++)
	lambda[count] = sqrt(lambda[count]);

      double lambdaTotal = lambda.norm1();
      double lambdaSum = 0.0;
      R_ = 0;
      for(int i = 0; i<lambda.dim(); i++)
	{
	  lambdaSum += lambda[i];
	  if(lambdaSum/lambdaTotal >= tol_)
	    {
	      // R is the number of modes to keep
	      R_ = i + 1;
	      SUNDANCE_ROOT_MSG2(verbosity_, "Number of modes kept: " + Teuchos::toString(R_));
	      break;
	    }
	}

      // Looking at PlayaSVD.cpp, Psi is a DenseSerialMatrix
      Playa::Vector<double> ej = Psi.domain().createMember();
      Playa::Vector<double> psiCoeff = Psi.range().createMember(); // These are the coefficient vectors
      psi_.resize(R_); // These are the pressure POD basis functions
      SUNDANCE_ROOT_MSG2(verbosity_, "Size of psi: " + Teuchos::toString(psi_.size()));

      // Get the Expr psi_r(x)
      CellFilter interior = new MaximalCellFilter();
      QuadratureFamily quad4 = new GaussianQuadrature(4);
      for(int r = 0; r<R_; r++)
	{
	  ej.zero();
	  ej[r] = 1.0;
	  psiCoeff.zero();
	  Psi.apply(ej,psiCoeff);
	  psi_[r] = new DiscreteFunction(ds_, serialToEpetra(psiCoeff)); //DiscreteFunction requires Epetra vectors
	  TEUCHOS_TEST_FOR_EXCEPTION( fabs(L2Norm(ds_.mesh(), interior, psi_[r], quad4) - 1.0) >= 1.0e-6,
				     runtime_error, "||psi["+Teuchos::toString(r)+"]|| = " + Teuchos::toString(L2Norm(ds_.mesh(), interior, psi_[r], quad4)) + " != 1");
	}

}


void PMB::generate_beta()
{
  // Set-up beta
  beta_.resize(M_);
  // Based off the value for R_, create an appropriate VectorSpace<double>
  VectorType<double> R_vecType = new SerialVectorType();
  VectorSpace<double> R_vecSpace = R_vecType.createEvenlyPartitionedSpace(MPIComm::self(), R_);
  for(int count = 0; count < M_; count++)
    beta_[count] = R_vecSpace.createMember();
  
  // c will hold c(t_m)
  Vector<double> c = R_vecSpace.createMember();
  Expr grad = gradient(ds_.mesh().spatialDim());

  LinearOperator<double> A(rcp(new DenseSerialMatrix(R_vecSpace, R_vecSpace)));
  RCP<DenseSerialMatrix> APtr = DenseSerialMatrix::getConcretePtr(A);

  for(int i = 0; i < R_; i++)
    {
      for(int j = 0; j < R_; j++)
	{
	  APtr->setElement(i,j,gradIP(psi_[i],psi_[j]));
	}
    }

  Expr nHat = CellNormalExpr(ds_.mesh().spatialDim(), "nHat");
  CellFilter boundary = new BoundaryCellFilter();
  QuadratureFamily quad = new GaussianQuadrature(6);

  
  DenseLUSolver solver;
  for(int time = 0; time < M_; time++)
    {
      t_.setParameterValue( time*deltat_ );
      c.zero();
      //      std::cout << "Here is q(t_m) " << forceTerm_ << std::endl;

      // Build c(t_m)
      for(int i = 0; i < R_; i++)
	{
	  //Expr integrand_boundary = nHat*(psi_[i]*(grad*pbar_)) + nHat*( psi_[i] * ((uRO_[time]*grad)*uRO_[time]) );
	  Expr integrand_boundary = nHat*( psi_[i] * ((uRO_[time]*grad)*uRO_[time]) );
	  FunctionalEvaluator integral_boundary = FunctionalEvaluator(ds_.mesh(), Integral(boundary, integrand_boundary, quad) );
	  
	  c[i] = -1.0*L2IP(grad*pbar_, grad*psi_[i]) - L2IP( (uRO_[time]*grad)*uRO_[time] , grad*psi_[i]) - L2IP(div(forceTerm_), psi_[i]) + integral_boundary.evaluate();
	}
      SUNDANCE_ROOT_MSG1(verbosity_, "Solving for time step " + Teuchos::toString(time) + " of " + Teuchos::toString(M_-1));
      SolverState<double> state = solver.solve(A,c,beta_[time]);
      TEUCHOS_TEST_FOR_EXCEPTION(state.finalState() != SolveConverged,
				 runtime_error, "solve failed");
    }
}


double PMB::gradIP(Expr f, Expr g)
  {
    CellFilter interior = new MaximalCellFilter();
    CellFilter boundary = new BoundaryCellFilter();
    QuadratureFamily quad = new GaussianQuadrature(6);
	  
    // mesh.spatialDim() returns n for nD
    int dim = ds_.mesh().spatialDim();

    // Define our differential operators; note Derivative(x=0)
    Expr grad = gradient(dim);
    Expr nHat = CellNormalExpr(dim,"nHat");

    Expr integrand_interior = (grad*f)*(grad*g);
    //Expr integrand_boundary = nHat*(f*(grad*g)); // Increased error from .04 to .11
    //Expr integrand_boundary = f*(nHat*(grad*g)); // Same result
    Expr integrand_boundary = f*0.0;
    
    FunctionalEvaluator IP = FunctionalEvaluator(ds_.mesh(), Integral(interior, integrand_interior, quad) - Integral(boundary, integrand_boundary, quad) ); 
    return (IP.evaluate());
  }


double PMB::L2IP(Expr f, Expr g)
  {
    CellFilter interior = new MaximalCellFilter();
    QuadratureFamily quad = new GaussianQuadrature(4);
	  
    FunctionalEvaluator IP = FunctionalEvaluator(ds_.mesh(), Integral(interior, f*g, quad));
    return (IP.evaluate());
  }

Array<Expr> PMB::get_pRO()
{

  Array<Expr> pRO(M_, pbar_);
  for(int time = 0; time < M_; time++)
    {
      for(int r = 0; r < R_; r++)
	{
	  pRO[time] = pRO[time] + beta_[time][r]*psi_[r];
	}
    }
  return pRO;
}
