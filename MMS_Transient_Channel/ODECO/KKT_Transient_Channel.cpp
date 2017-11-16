#include "KKT_Transient_Channel.hpp"


KKT_Transient_Channel::KKT_Transient_Channel(string POD_DataDir, Mesh spatialMesh, Mesh timeMesh, sensorData dataClass, double tFinal, int nSteps, int quadOrder, int verbosity)
  : KKTBase(POD_DataDir, spatialMesh, timeMesh, tFinal, nSteps, verbosity),
    dataClass_(dataClass),
    quad_(new GaussianQuadrature(quadOrder))
{}

void KKT_Transient_Channel::stateEqn()
{
  CellFilter verts = new DimensionalCellFilter(0);
  CellFilter left = verts.coordSubset(0,0.0); // hard coded to start at time = 0
  // H(T,alpha) is implemented in sundance as T*alpha*alpha
  stateEqn_ = Integral(interior_, lambdaHat_*(dt_*alpha_ - A_*alpha_ - T_*alpha_*alpha_ -b_), quad_);
  stateBC_ = EssentialBC(left, lambdaHat_*(alpha_-p_), quad_);

}

void KKT_Transient_Channel::adjointEqn()
{

  integralOperator S = dataClass_.get_S();
  Expr vstar = dataClass_.get_vstar();
  Array<Point> sensorLocations = dataClass_.get_locations();
  int Ns = vstar.size();  

  // Build Su, where S_i(u) = S_i(uB) + sum( alpha[r] S_i( phi[r] ) )
  Expr Su;
  for(int i = 0; i < Ns; i++)
    {
      // Calculate S_i(uB)
      CellFilter location = getPoint(sensorLocations[i]);
      Expr Si = S.staticDetect(location, uB_);
      for(int r = 0; r < Ru_; r++)
  	{
  	  Si = Si + alpha_[r]*S.staticDetect(location,phi_[r]);
  	}
      Su.append(Si);
    }

  // The purpose of this section of code is to calculate (Su - vstar, S(PHI*alphaHat) )
  // For the moment let the IP be the dot product
  Expr SPHIaHat;
  for(int i = 0; i < Ns; i++)
    {
      CellFilter location = getPoint(sensorLocations[i]);
      Expr Si = 0.0;
      for(int r = 0; r < Ru_; r++)
  	{
  	  Si = Si + alphaHat_[r]*S.staticDetect(location,phi_[r]);
  	}
      SPHIaHat.append(Si);
    }
      
  adjointEqn_ = Integral(interior_, (Su - vstar)*SPHIaHat - ( alphaHat_*(dt_*lambda_ - At_*lambda_) - lambda_*(T_*alphaHat_*alpha_ + T_*alpha_*alphaHat_) ), quad_);

  CellFilter verts = new DimensionalCellFilter(0);
  CellFilter right = verts.coordSubset(0,tFinal_);
  adjointBC_ = EssentialBC(right, alphaHat_*lambda_, quad_);
}

void KKT_Transient_Channel::designEqn(double eta)
{
  /* design equation and BC */
  /* -- the (alpha')*(alphaHat') term enforces constancy of alpha in time */
  CellFilter verts = new DimensionalCellFilter(0);
  CellFilter left = verts.coordSubset(0,0.0); // hard coded to start at time = 0
  designEqn_ = Integral(interior_, eta*(dt_*p_)*(dt_*pHat_), quad_);
  designBC_ = EssentialBC(left, pHat_*lambda_, quad_);
}



void KKT_Transient_Channel::regularizationTerm(double eta)
{
  Expr regTerm = 0.0;
  // summation over i
  for(int i = 0; i < Ru_; i++)
    {
      FunctionalEvaluator temp(spatialMesh_, Integral(interior_, uB_*phi_[i], quad_));
      regTerm = regTerm + alphaHat_[i]*temp.evaluate();
    }

  // // summation over j and k
  // for(int j = 0; j < Ru; j++)
  // 	{
  // 	  for(int k = 0; k < Ru; k++)
  // 	    {
  // 	      FunctionalEvaluator temp(spatialMesh, Integral(interior, phi[j]*phi[k],quad));
  // 	      cout << "(phi[" << j << "], phi[" << k << "]) = " << temp.evaluate() << endl;
  // 	      regTerm = regTerm + alpha[j]*alphaHat[k]*temp.evaluate();
  // 	    }
  // 	}

  // Summation over j and k simplifying due to orthogonality of the basis functions
  for(int j = 0; j < Ru_; j++)
    {
      regTerm = regTerm + alpha_[j]*alphaHat_[j];
    }

  // Multiply by the constant term eta
  regTerm = regTerm*eta;
  adjointEqn_ = adjointEqn_ + Integral(interior_, regTerm, quad_);
}






Expr KKT_Transient_Channel::solve(string solverXML, double eta_design, double eta_reg)
{
  Array<BasisFamily> ODECO_basisArray(3*Ru_);
  for (int i=0; i<ODECO_basisArray.size(); i++)
    ODECO_basisArray[i]=time_basis_;
  DiscreteSpace ODECO_DS(timeMesh_, ODECO_basisArray, epetraVecType_);
  // NonlinearSolver will override the intial guess with the final solution
  Expr U0 = new DiscreteFunction(ODECO_DS, 0.0);

  // Due to the fact that the initial state of the approximation is terrible for high Re,
  // I will give the first Ru_ components of U0 the underlying Vector<double> for alpha
  // These values will be obtained via projeciton
  // This seemed to have no effect
  /*
  Expr u0 = List(1.0,0.0);
  VectorType<double> Ru_vecType = new SerialVectorType();
  VectorSpace<double> Ru_vecSpace = Ru_vecType.createEvenlyPartitionedSpace(MPIComm::self(),Ru_);
  Vector<double> a0 = Ru_vecSpace.createMember();
  for(int i = 0; i < Ru_; i++)
    {
      FunctionalEvaluator IP = FunctionalEvaluator(spatialMesh_, Integral(interior_,
									  (u0 - uB_)*phi_[i],
									  quad_));
      a0[i] = IP.evaluate();
    }
  Vector<double> U0_Vec = getDiscreteFunctionVector(U0);
  std::cout << "Size of U0_Vec: " << U0_Vec.dim() << std::endl;
  // Need to give the first Ru_ components of a0's vec to U0_Vec
  for(int r = 0; r < Ru_; r++)
    U0_Vec[r] = a0[r];
  */

  // Next thought: Take the alpha from the inverse problem and use that to build U0_Vec
  Vector<double> U0_Vec = getDiscreteFunctionVector(U0);
  string filename = "/home/sirush/PhDResearch/MMS_Transient_Channel/InverseProblem/Results/SingleParameterSpace/tol0.999/Re128/nx25nt25/alphaROM.txt";
  std::ifstream alpha_reader(filename, std::ios::in);
  TEUCHOS_TEST_FOR_EXCEPTION(!alpha_reader,
			     std::runtime_error,
			     "could not open file " << filename << " for reading");
  int basisOrder;
  alpha_reader >> basisOrder;
  int temp;
  alpha_reader >> temp; // Ru is the number of elements in vector-valued function b
  int numOfVectors; // equal to nSteps+1
  alpha_reader >> numOfVectors;
  string alias;
  alpha_reader >> alias;

  // Check to make sure that the time-dependent only functions are using Lagrange(1)
  TEUCHOS_TEST_FOR_EXCEPTION(basisOrder != 1,
			     std::runtime_error,
			     "The order of the basis functions for time-dependent only functions is not Lagrange(1)");


  // Create the time DiscreteSpace
  Array<BasisFamily> time_basisArray_(Ru_);
  for(int i = 0; i < Ru_; i++)
    time_basisArray_[i] = time_basis_;
  //  time_DS_ = new DiscreteSpace(timeMesh_, time_basisArray, epetraVecType_);
  DiscreteSpace time_DS_(timeMesh_, time_basisArray_, epetraVecType_);
  
  // Establish the values for U0
  for(int tIndex=0; tIndex < nSteps_+1; tIndex++)
    {
      for(int r = 0; r < Ru_; r++)
	// I think that the problem is that there is an implied third runner
	// for which function (alpha, lambda, p)  is being referenced.
	// Thus I propose that the runners go t, k, r instead of t, r, k
	alpha_reader >> U0_Vec[r + 3*Ru_*tIndex];
    }  

  stateEqn();
  adjointEqn();
  regularizationTerm(eta_reg);
  designEqn(eta_design);


  // Combine the equations and BCs to form the KKT system 
  Expr eqn = stateEqn_ + adjointEqn_ + designEqn_;
  Expr bc = stateBC_ + adjointBC_ + designBC_;

  // create the NLP
  NonlinearProblem NLP(timeMesh_, eqn, bc, List(lambdaHat_, alphaHat_, pHat_),
		       List(alpha_, lambda_, p_), U0, epetraVecType_);
  
  NonlinearSolver<double> solver 
    = NonlinearSolverBuilder::createSolver(solverXML);

  // Solve
  cout << "Solving problem..." << endl;
  SolverState<double> state = NLP.solve(solver);    

  // Verify the solver converged
  TEUCHOS_TEST_FOR_EXCEPTION(state.finalState() != SolveConverged,
			     std::runtime_error,
			     "Nonlinear solve failed to converge: message="
			     << state.finalMsg());

  // This is what I had originally; I think it is just storing
  // copies of U0, not elements of it since it is an Expr
  // alphaOPT_ = U0[0];
  // for(int r=1; r < Ru_; r++)
  //   alphaOPT_.append(U0[r]);

  // This did not work; says alphaEx has 2 elements and alphaOPT has 6
  //alphaOPT_ = U0;

  // Attempting to create alphaOPT_ from the first Ru_*(nSteps_+1) elements
  // getDiscreteFunctionVector(U0).dim()
  U0_Vec = getDiscreteFunctionVector(U0);
  Array<BasisFamily> time_basisArray(Ru_);
  for(int i = 0; i < Ru_; i++)
    time_basisArray[i] = time_basis_;
  DiscreteSpace time_DS(timeMesh_, time_basisArray, epetraVecType_);
  alphaOPT_ = new DiscreteFunction(time_DS, 0.0, "alphaOPT_");
  // The get() links to the DiscreteFunction's actual vector (shallow)
  Vector<double> alphaOPT_Vec = getDiscreteFunctionVector(alphaOPT_);
  for(int tIndex=0; tIndex < nSteps_+1; tIndex++)
    {
      for(int r = 0; r < Ru_; r++)
	// KL says this ordering should work, but it's not
	//alphaOPT_Vec[r+tIndex*Ru_] = U0_Vec[3*r+3*Ru_*tIndex];

	// I think that the problem is that there is an implied third runner
	// for which function (alpha, lambda, p)  is being referenced.
	// Thus I propose that the runners go t, k, r instead of t, r, k
	alphaOPT_Vec[r+tIndex*Ru_] = U0_Vec[r + 3*Ru_*tIndex];
    }

  // Array<Expr> temp(nSteps_+1,uB_);
  // uOPT_ = temp;
  for(int time = 0; time < nSteps_+1; time++)
    {
      for(int r = 0; r < Ru_; r++)
  	{
  	  uOPT_[time] = uOPT_[time] + alphaOPT_Vec[r + time*Ru_]*phi_[r];
  	}
    }

  // This is correct, but results in alphaOPT_ having the same underlying
  // Vector for each component
  // alphaOPT_ = U0[0];
  // for(int i=1; i< Ru_; i++)
  //   {
  //     // From the file 2DAnalytical.... in the KLWay folder
  //     alphaOPT_.append(U0[i]);
  //   }

  // Each component of U0 shares a single underlying vector of size 3*(Ru_*(nSteps+1))
  // for(int i=0; i< Ru_*3; i++)
  //   {
  //     std::cout << "This is the size of  U0["<<i<<"]'s vector: "
  // 		<< getDiscreteFunctionVector(U0[i]).dim() << std::endl;
  //   }

  return U0;
}

Array<double> KKT_Transient_Channel::errorCheck(Expr uExact, Expr t)
{
  // Based off the value for Ru, create an appropriate VectorSpace<double>
  VectorType<double> Ru_vecType = new SerialVectorType();
  VectorSpace<double> Ru_vecSpace = Ru_vecType.createEvenlyPartitionedSpace(MPIComm::self(), Ru_);

  // Find the Array<Vector> representing alphaExact(t_m) for m=0:nSteps
  Array<Vector<double> > ipVec(nSteps_+1);
  for(int count = 0; count<ipVec.length(); count++)
    ipVec[count] = Ru_vecSpace.createMember();

  //Needed for the integral
  //QuadratureFamily quad = new GaussianQuadrature(2);

  //Find the vectors with the formula alphaEx_r(t_m) = (uExact(t_m,x,y) - uB, phi_[r] )
  double tInit = 0.0;
  double deltat = (tFinal_ - tInit)/nSteps_;
  for(int tIndex=0; tIndex < ipVec.length(); tIndex++)
    {
      t.setParameterValue(tInit+tIndex*deltat);
      for(int r=0; r<Ru_; r++)
	{
	  FunctionalEvaluator ExactEvaluator = FunctionalEvaluator(spatialMesh_, Integral(interior_, (uExact-uB_)*phi_[r], quad_));
	  ipVec[tIndex][r] = ExactEvaluator.evaluate();
	}
    }

  // Turn alphaExVec into a DiscreteFunction based off the values in ipVec
  // Create the time DiscreteSpace
  Array<BasisFamily> time_basisArray(Ru_);
  for(int i = 0; i < Ru_; i++)
    time_basisArray[i] = time_basis_;
  DiscreteSpace time_DS(timeMesh_, time_basisArray, epetraVecType_);
  Expr alphaExact = new DiscreteFunction(time_DS, 0.0, "alphaExact");
  // The get() links to the DiscreteFunction's actual vector (shallow)
  Vector<double> alphaExact_Vec = getDiscreteFunctionVector(alphaExact);
  for(int tIndex=0; tIndex < nSteps_+1; tIndex++)
    {
      for(int r = 0; r < Ru_; r++)
	alphaExact_Vec[r+tIndex*Ru_] = ipVec[tIndex][r];
    }

  // Compare alphaExact and alphaOPT
  double alphaAbsError = L2Norm(timeMesh_, interior_, alphaExact-alphaOPT_, quad_);
  double alphaRelError = alphaAbsError/L2Norm(timeMesh_, interior_, alphaExact, quad_);


  /*
   * Begin calculating the error comparison for uOPT against uExact
   * Since we cannot have both time and x be CoordExpr(0), we will discretize
   * in time. This will require getting the underlying Vector<double>
   * from alphaOPT_, which will have Ru_*(nSteps_+1) elements
   */
  SUNDANCE_ROOT_MSG2(verbosity_, "(Ru_ = " << Ru_ << ")*(nSteps_+1="
		     << nSteps_+1 << ") = " << Ru_*(nSteps_+1) );
  TEUCHOS_TEST_FOR_EXCEPTION(getDiscreteFunctionVector(alphaOPT_).dim() != alphaExact_Vec.dim(),
			     runtime_error,
			     "The size of the underlying vector for alphaOPT_ ("
			     << getDiscreteFunctionVector(alphaOPT_).dim()
			     << " ) is not equal to the size of the underlying vector for "
			     << " alphaExact (" << alphaExact_Vec.dim() << ")");
  SUNDANCE_ROOT_MSG1(verbosity_, "Size of ||alphaEx|| "
		     << L2Norm(timeMesh_, interior_, alphaExact, quad_) );

  Vector<double> alphaOPT_Vec = getDiscreteFunctionVector(alphaOPT_);



  // Array<Expr> uOPT_(nSteps_+1, uB_);
  // for(int time = 0; time < nSteps_+1; time++)
  //   {
  //     for(int r = 0; r < Ru_; r++)
  // 	{
  // 	  uOPT_[time] = uOPT_[time] + alphaOPT_Vec[r + time*Ru_]*phi_[r];
  // 	}
  //   }

  VectorType<double> timeVecType = new SerialVectorType();
  VectorSpace<double> timeVecSpace = timeVecType.createEvenlyPartitionedSpace(MPIComm::self(), nSteps_+1.0);

  //Needed for the integral
  QuadratureFamily quadU = new GaussianQuadrature(6);

  Vector<double> velocityAbsError = timeVecSpace.createMember();
  Vector<double> velocityRelError = timeVecSpace.createMember();
  for(int time = 0; time < nSteps_+1; time++)
    {
      t.setParameterValue(tInit+time*deltat);
      velocityAbsError[time] = L2Norm(spatialMesh_, interior_, uExact - uOPT_[time], quadU);
      velocityRelError[time] = velocityAbsError[time]/L2Norm(spatialMesh_, interior_, uExact, quadU);
      //cout << "||uExact(t="<<time<<",x,y)||_2 = " << L2Norm(spatialMesh_, interior_, uExact, quadU) << endl;
    }

  // Convert the estimations of the error for velocity into single double values
  // by taking the two norm of each
  double aggregateVelAbsError = velocityAbsError.norm2();
  double aggregateVelRelError = velocityRelError.norm2();
  
  return tuple(alphaAbsError, alphaRelError, aggregateVelAbsError, aggregateVelRelError);
  
}
