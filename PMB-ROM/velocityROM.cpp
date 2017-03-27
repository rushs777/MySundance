#include "velocityROM.hpp"

velocityROM::velocityROM(const string snapshotFilename, string nlParamFile, const DiscreteSpace& ds, Expr u0, Expr forceTerm, Expr t, int nSteps, double deltat, double tolerance, int verbosity)
  : snapshotFilename_(snapshotFilename),
    nlParamFile_(nlParamFile),
    ds_(ds),
    u0_(u0),
    forceTerm_(forceTerm),
    t_(t),
    nSteps_(nSteps), deltat_(deltat),
    tol_(tolerance),
    verbosity_(verbosity)
{
  //  t_ = new Sundance::Parameter(0.0);
  t_.setParameterValue(0.0);
  Expr x = new CoordExpr(0,"x");
  Expr y = new CoordExpr(1,"y");
  /*
  double nu = 1.0;
  forceTerm_ = List((-36*cos(y)*sin(t_)*sin(x) + (((44 + 30*cos(t_) + 8*cos(4*t_) + 30*cos(5*t_))*
						   cos(x) + 9*((3 + 2*cos(t_) + 2*cos(5*t_))*cos(3*x) + cos(5*x)))*
						  pow(sin(x),3) + 9*cos(6*t_)*cos(2*x)*pow(sin(2*x),3))*pow(sin(y),2) - 
	       3*(8*nu*cos(2*t_)*(-1 + 2*cos(2*x)) + 6*nu*cos(3*t_)*(-1 + 5*cos(4*x)) + 
		  8*sin(2*t_)*pow(sin(x),2) + 9*sin(3*t_)*pow(sin(2*x),2))*sin(2*y))/36.,
	      (6*nu*(2*cos(2*t_)*(-1 + 2*cos(2*y))*sin(2*x) + 
		     3*cos(3*t_)*(-4 + 5*cos(2*y))*sin(4*x)) - 18*cos(x)*sin(t_)*sin(y) + 
	       3*(4*sin(2*t_)*sin(2*x) + 9*sin(3*t_)*sin(4*x))*pow(sin(y),2) + 
	       2*cos(y)*(cos(2*t_)*(4*cos(2*t_) - 3*cos(3*t_)*(-3 - 6*cos(2*x) + cos(4*x)))*
			 pow(sin(x),2) + 9*pow(cos(3*t_),2)*pow(sin(2*x),2))*pow(sin(y),3))/
	      18.);
  */
}

void velocityROM::initialize()
{
  Wprime_ = snapshotToMatrix(snapshotFilename_, nSteps_, ds_.mesh());
  SUNDANCE_ROOT_MSG2(verbosity_, "Size of W: " << Wprime_.range().dim() << " by " << Wprime_.domain().dim());

  // Calculate ubar(x)
  Vector<double> ubarVec = Wprime_.range().createMember();
  Vector<double> ones = Wprime_.domain().createMember();
  //ones.setToConstant(0.0);
  ones.setToConstant(1.0);
  Wprime_.apply(ones, ubarVec);
  ubarVec *= (1.0/ (nSteps_+1.0) );
  ubar_ = new DiscreteFunction(ds_, serialToEpetra(ubarVec));

  // Make W be W'
  RCP<DenseSerialMatrix> WPtr = DenseSerialMatrix::getConcretePtr(Wprime_);
  double Wij;
  for(int j = 0; j < Wprime_.domain().dim(); j++)
    {
      for(int i = 0; i<Wprime_.range().dim(); i++)
	{
	  Wij = WPtr->getElement(i,j);
	  WPtr->setElement(i,j,Wij-ubarVec[i]);
	}  
    }
      // Perform the POD of the matrix W'
      Playa::LinearOperator<double> U;
      Playa::LinearOperator<double> Phi;
      Playa::Vector<double> sigma;


      // W' and mesh need to be defined
      SUNDANCE_ROOT_MSG1(verbosity_, "Entering POD for velocity");
      POD(Wprime_,sigma,U,Phi,ds_,verbosity_);
      SUNDANCE_ROOT_MSG2(verbosity_, "POD finished for velocity");

      Vector<double> lambda = sigma.copy();
      for(int count = 0; count < sigma.dim(); count++)
	lambda[count] = sqrt(lambda[count]);

      double lambdaTotal = lambda.norm1();
      double lambdaSum = 0.0;
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

      // Looking at PlayaSVD.cpp, Phi is a DenseSerialMatrix
      Playa::Vector<double> ej = Phi.domain().createMember();
      Playa::Vector<double> phiCoeff = Phi.range().createMember(); // These are the coefficient vectors
      phi_.resize(R_); // These are the velocity POD basis functions
      SUNDANCE_ROOT_MSG2(verbosity_, "Size of phi: " + Teuchos::toString(phi_.size()));

      // Get the Expr phi_r(x)
      CellFilter interior = new MaximalCellFilter();
      QuadratureFamily quad4 = new GaussianQuadrature(6);
      for(int r = 0; r<R_; r++)
	{
	  ej.zero();
	  ej[r] = 1.0;
	  phiCoeff.zero();
	  Phi.apply(ej,phiCoeff);
	  phi_[r] = new DiscreteFunction(ds_, serialToEpetra(phiCoeff)); //DiscreteFunction requires Epetra vectors
	  TEUCHOS_TEST_FOR_EXCEPTION( fabs(L2Norm(ds_.mesh(), interior, phi_[r], quad4) - 1.0) >= 1.0e-6,
				     runtime_error, "||phi["+Teuchos::toString(r)+"]|| = " + Teuchos::toString(L2Norm(ds_.mesh(), interior, phi_[r], quad4)) + " != 1");
	}

      // Generate alpha[0] based off of u0_ and ubar_
      // (u0_ - ubar_, phi[i]) = alpha[0][i]
      alpha_.resize(nSteps_+1);

      VectorType<double> R_vecType = new SerialVectorType();
      VectorSpace<double> R_vecSpace = R_vecType.createEvenlyPartitionedSpace(MPIComm::self(), R_);
      alpha_[0] = R_vecSpace.createMember();
      
      for(int i = 0; i < R_; i++)
	{
	  FunctionalEvaluator IP = FunctionalEvaluator(ds_.mesh(), Integral(interior, (u0_ - ubar_)*phi_[i], quad4));
	  //FunctionalEvaluator IP = FunctionalEvaluator(ds_.mesh(), Integral(interior, u0_*phi_[i], quad4));
	  alpha_[0][i] = IP.evaluate();
	}

}


void velocityROM::generate_alpha()
{
  SUNDANCE_ROOT_MSG1(verbosity_, "Creating MMSQuadODE");
  // Create the nonlinear operator for solving our nonlinear ODE
  MMSQuadODE f(phi_, ubar_, forceTerm_, t_, deltat_, ds_.mesh(), false, verbosity_);
  f.initialize();

  SUNDANCE_ROOT_MSG1(verbosity_, "Creating NLO");
  MyNLO* prob = new MyNLO(f, deltat_);
  NonlinearOperator<double> F = prob;
  if(verbosity_>=3)
    F.setVerb(verbosity_);


  SUNDANCE_ROOT_MSG1(verbosity_, "Creating solver");
  // create the Newton-Armijo solver
  ParameterXMLFileReader reader(nlParamFile_);
  ParameterList params = reader.getParameters();
  const ParameterList& solverParams = params.sublist("NewtonArmijoSolver");
  
  //Next line possible since DenseLUSolver and LinearSolver both inherit from LinearSolverBase
  LinearSolver<double> linearSolver(rcp(new DenseLUSolver())); 
  NewtonArmijoSolver<double> nonlinearSolver(solverParams, linearSolver);


  // Establish alpha_(t_init)
  SUNDANCE_ROOT_MSG1(verbosity_, "Running nonlinear solves...");
  prob->set_uPrev(alpha_[0]);
  for(int time = 1; time < alpha_.length(); time++)
    {
      SUNDANCE_ROOT_MSG1(verbosity_, "Nonlinear Solve at time step " + Teuchos::toString(time) + " of " + Teuchos::toString(nSteps_));
      prob->set_tPrev( (time-1.0)*deltat_ );
      SolverState<double> state = nonlinearSolver.solve(F, alpha_[time]);
      
      if(verbosity_>=3)
	cout << "alpha[" << Teuchos::toString(time) << "]: " <<  alpha_[time];
      
      TEUCHOS_TEST_FOR_EXCEPTION(state.finalState() != SolveConverged,
				 runtime_error, "solve failed");
      prob->set_uPrev(alpha_[time]);
    }

}


Array<Expr> velocityROM::get_uRO()
{
  Array<Expr> uRO(nSteps_+1, ubar_);
  for(int time = 0; time < nSteps_+1; time++)
    {
      for(int r = 0; r < R_; r++)
	{
	  uRO[time] = uRO[time] + alpha_[time][r]*phi_[r];
	}
    }
  return uRO;
}
