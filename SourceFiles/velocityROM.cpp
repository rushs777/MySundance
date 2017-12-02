#include "velocityROM.hpp"

velocityROM::velocityROM(const string POD_base_dir, string nlParamFile,
			 const DiscreteSpace& ds, Expr u0, Expr forceTerm, Expr t,
			 int nSteps, double deltat, double nu, double tolerance, int verbosity)
  : POD_base_dir_(POD_base_dir),
    nlParamFile_(nlParamFile),
    ds_(ds),
    u0_(u0),
    forceTerm_(forceTerm),
    t_(t),
    nSteps_(nSteps),
    deltat_(deltat),
    nu_(nu),
    tol_(tolerance),
    verbosity_(verbosity)
{
  t_.setParameterValue(0.0);
}

void velocityROM::initialize(const string snapshotFilePrefix)
{
  // At this point, the matrix is really W
  LinearOperator<double> W = snapshotToMatrix(snapshotFilePrefix, nSteps_, ds_.mesh());
  SUNDANCE_ROOT_MSG2(verbosity_, "Size of W: " << W.range().dim() << " by " << W.domain().dim());

  // calculate uB_(x)
  uB_ = timeMeanFunctionGenerator(W,ds_);


  // Read in the POD from file
  SUNDANCE_ROOT_MSG1(verbosity_, "Reading in the reduced-order basis for velocity");
  string POD_basis_fileprefix = POD_base_dir_ + "/POD_basis";

  R_ = 10000;
  for(int r = 0; r < R_; r++)
    {
      try
	{
	  phi_.push_back( readSnap(POD_basis_fileprefix, r, ds_.mesh() ) );
	}
      catch (std::runtime_error& e)
	{
	  R_ = r;
	}
    }
  
  SUNDANCE_ROOT_MSG1(verbosity_, "Found " << R_ << " reduced-order basis functions for the given resoultion"); 


  // Generate alpha[0] based off of u0_ and uB_
  // (u0_ - uB_, phi[i]) = alpha[0][i]
  CellFilter interior = new MaximalCellFilter();
  QuadratureFamily quad = new GaussianQuadrature(6);
  
  alpha_.resize(nSteps_+1);

  VectorType<double> R_vecType = new SerialVectorType();
  VectorSpace<double> R_vecSpace = R_vecType.createEvenlyPartitionedSpace(MPIComm::self(), R_);
  alpha_[0] = R_vecSpace.createMember();
      
  for(int i = 0; i < R_; i++)
    {
      FunctionalEvaluator IP = FunctionalEvaluator(ds_.mesh(), Integral(interior, (u0_ - uB_)*phi_[i], quad));
      alpha_[0][i] = IP.evaluate();
    }

}

// Resume by getting uRO working for this new initialization
void velocityROM::generate_alpha()
{
  SUNDANCE_ROOT_MSG1(verbosity_, "Creating NSEProjectedODE");
  // Create the nonlinear operator for solving our nonlinear ODE
  //NSEProjectedODE f(phi_, ubar_, forceTerm_, t_, deltat_, ds_.mesh(), false, verbosity_);
  RCP<NSEProjectedODE> f = rcp(new NSEProjectedODE(phi_, uB_, forceTerm_, t_, deltat_, nu_,
						   ds_.mesh(), false, verbosity_));
  f->initialize();

  SUNDANCE_ROOT_MSG1(verbosity_, "Creating NLO");
  MyNLO* prob = new MyNLO(f, deltat_);
  NonlinearOperator<double> F = prob;
  if(verbosity_>=4)
    F.setVerb(verbosity_);


  SUNDANCE_ROOT_MSG1(verbosity_, "Creating solver.......");
  // create the Newton-Armijo solver
  ParameterXMLFileReader reader(nlParamFile_);
  ParameterList params = reader.getParameters();
  const ParameterList& solverParams = params.sublist("NewtonArmijoSolver");
  
  //Next line possible since DenseLUSolver and LinearSolver both inherit from LinearSolverBase
  LinearSolver<double> linearSolver(rcp(new DenseLUSolver())); 
  NewtonArmijoSolver<double> nonlinearSolver(solverParams, linearSolver);


  // Establish alpha_(t_init)
  SUNDANCE_ROOT_MSG1(verbosity_, "Running nonlinear solves.......");
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

void velocityROM::write_b(const string filename)
{
  RCP<NSEProjectedODE> f = rcp(new NSEProjectedODE(phi_, uB_, forceTerm_, t_, deltat_,
						   nu_, ds_.mesh(), false, verbosity_));
  f->initialize();

  // // Assumes that tInit = 0.0
  // double tInit = 0.0;
  // std::ofstream writer(filename);
  // writer << "1" << endl;
  // writer << R_ << endl;
  // writer << nSteps_ + 1 << endl;
  // writer << "b" << endl;
  // Vector<double> b_ti = alpha_[0].copy();

  // for(int i = 0; i < nSteps_+1; i++)
  //   {
      
  //     b_ti = f->evalForceTerm(tInit + i*deltat_);
  //     for(int r = 0; r < R_; r++)
  // 	{
  // 	  writer << b_ti[r] << endl;
  // 	}
  //   }

  f->write_b(nSteps_);
  
}

void velocityROM::write_alpha(const string filename)
{
  // Assumes that tInit = 0.0
  double tInit = 0.0;
  std::ofstream writer(filename);
  int lagrangeBasis = 1;
  writer << lagrangeBasis << endl;
  writer << R_ << endl; // Number of components for each time step
  writer << nSteps_ + 1 << endl;
  string alias = "alphaROM";
  writer << alias << endl;

  for(int timeStep = 0; timeStep < nSteps_ + 1; timeStep++)
    {
      for(int r = 0; r < R_; r++)
	{
	  writer << alpha_[timeStep][r] << endl;
	}
    }
}


Array<Expr> velocityROM::get_uRO()
{
  Array<Expr> uRO(nSteps_+1, uB_);
  for(int time = 0; time < nSteps_+1; time++)
    {
      for(int r = 0; r < R_; r++)
	{
	  uRO[time] = uRO[time] + alpha_[time][r]*phi_[r];
	}
    }
  return uRO;
}


Vector<double> velocityROM::alphaErrorCheck(Expr uExact)
{
  // For the purposes of comparison, start calculating alpha "exactly"    
  VectorType<double> timeVecType = new SerialVectorType();
  VectorSpace<double> timeVecSpace = timeVecType.createEvenlyPartitionedSpace(MPIComm::self(), nSteps_+1.0);

  // Based off the value for Ru, create an appropriate VectorSpace<double>
  VectorType<double> R_vecType = new SerialVectorType();
  VectorSpace<double> R_vecSpace = R_vecType.createEvenlyPartitionedSpace(MPIComm::self(), R_);


  // Find the exact alphas
  Array<Vector<double> > alphaEx(nSteps_+1);
  for(int count = 0; count<alphaEx.length(); count++)
    alphaEx[count] = R_vecSpace.createMember();

  //Needed for the integral
  CellFilter interior = new MaximalCellFilter();
  QuadratureFamily quad = new GaussianQuadrature(6);


  //Find the alphaEx with the formula alphaEx_r(t_m) = (uExact(t_m,x,y), phi_[r] )
  double tInit = 0.0;
  for(int tIndex=0; tIndex < alphaEx.length(); tIndex++)
    {
      t_.setParameterValue(tInit+tIndex*deltat_);
      for(int r=0; r<R_; r++)
	{
	  FunctionalEvaluator ExactEvaluator = FunctionalEvaluator(ds_.mesh(), Integral(interior, (uExact-uB_)*phi_[r], quad));
	  alphaEx[tIndex][r] = ExactEvaluator.evaluate();
	}
    }

  // Note alphaROM is called alpha_
  Vector<double> alphaError = timeVecSpace.createMember();
  for(int i=0; i < alpha_.length(); i++)
    {
      alphaError[i] = (alphaEx[i] - alpha_[i]).norm2();
      if(verbosity_>=2)
	cout << "Error for alpha(t=" << i << "): " << alphaError[i] << endl;

      if(verbosity_>=3)
	{
	  cout << "Exact alpha(t=" << i << "): " << endl << alphaEx[i] << endl;	
	  cout << "Approximate alpha(t=" << i << "): " << endl << alpha_[i] << endl << endl;
	}
    }

  return alphaError;
  
}


Vector<double> velocityROM::velocityErrorCheck(Expr uExact)
{
  // For the purposes of comparison, start calculating alpha "exactly"    
  VectorType<double> timeVecType = new SerialVectorType();
  VectorSpace<double> timeVecSpace = timeVecType.createEvenlyPartitionedSpace(MPIComm::self(), nSteps_+1.0);

  //Needed for the integral
  CellFilter interior = new MaximalCellFilter();
  QuadratureFamily quad = new GaussianQuadrature(6);

  Vector<double> velocityError = timeVecSpace.createMember();
  Array<Expr> uRO = get_uRO();

  double tInit = 0.0;
  for(int time = 0; time < nSteps_+1; time++)
    {
      t_.setParameterValue(tInit+time*deltat_);
      velocityError[time] = L2Norm(ds_.mesh(), interior, uExact - uRO[time], quad);
      double uExactNorm = L2Norm(ds_.mesh(), interior, uExact, quad);
      /* print the relative error */
      SUNDANCE_ROOT_MSG2(verbosity_, "Relative Error for uRO at time " << t_.getParameterValue()
			 << " = " << velocityError[time]/uExactNorm);
    }

  return velocityError;
}
