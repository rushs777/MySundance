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

  return U0;
}

double KKT_Transient_Channel::errorCheck(Expr alphaOPT, Expr uExact, Expr t)
{
  // Based off the value for Ru, create an appropriate VectorSpace<double>
  VectorType<double> Ru_vecType = new SerialVectorType();
  VectorSpace<double> Ru_vecSpace = Ru_vecType.createEvenlyPartitionedSpace(MPIComm::self(), Ru_);

  // Find the Array<Vector> representing alphaExact(t_m) for m=0:nSteps
  Array<Vector<double> > ipVec(nSteps_+1);
  for(int count = 0; count<ipVec.length(); count++)
    ipVec[count] = Ru_vecSpace.createMember();

  //Needed for the integral
  QuadratureFamily quad = new GaussianQuadrature(6);

  //Find the vectors with the formula alphaEx_r(t_m) = (uExact(t_m,x,y) - uB, phi_[r] )
  double tInit = 0.0;
  double deltat = (tFinal_ - tInit)/nSteps_;
  for(int tIndex=0; tIndex < ipVec.length(); tIndex++)
    {
      t.setParameterValue(tInit+tIndex*deltat);
      for(int r=0; r<Ru_; r++)
	{
	  FunctionalEvaluator ExactEvaluator = FunctionalEvaluator(spatialMesh_, Integral(interior_, (uExact-uB_)*phi_[r], quad));
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

  // Compare alphaROM and alphaOPT
  double alphaError = L2Norm(timeMesh_, interior_, alphaExact-alphaOPT, quad);
  return alphaError;

  //alphaExact_Vec.space().numLocalElements()
  
  // VectorType<double> timeVecType = new SerialVectorType();
  // VectorSpace<double> timeVecSpace = timeVecType.createEvenlyPartitionedSpace(MPIComm::self(), nSteps_+1.0);
  // Vector<double> alphaError = timeVecSpace.createMember();  
  // for(int i=0; i < alphaOPT.length(); i++)
  //   {
  //     alphaError[i] = (alphaEx[i] - alphaOPT[i]).norm2();
  //     if(verbosity_>=2)
  // 	cout << "Error for alpha(t=" << i << "): " << alphaError[i] << endl;

  //     if(verbosity_>=3)
  // 	{
  // 	  cout << "Exact alpha(t=" << i << "): " << endl << alphaEx[i] << endl;	
  // 	  cout << "KKT alpha(t=" << i << "): " << endl << alphaOPT[i] << endl << endl;
  // 	}
  //   }



  // Legacy: double KKT_Transient_Channel::errorCheck(string ROM_base_dir, Expr alphaOPT)
  // Read in alphaROM(t) from the ROM code 
  // string alphaROM_filename = ROM_base_dir + "alphaROM.txt";
  // std::ifstream alphaROM_reader(alphaROM_filename, std::ios::in);
  // TEUCHOS_TEST_FOR_EXCEPTION(!alphaROM_reader, std::runtime_error,
  // 			     "could not open file"
  // 			     << alphaROM_filename << " for reading");

  // int basisOrder;
  // alphaROM_reader >> basisOrder;
  // int numOfElements;
  // alphaROM_reader >> numOfElements; // This is Ru
  // int numOfVectors;
  // alphaROM_reader >> numOfVectors;
  // string alias;
  // alphaROM_reader >> alias;


  // TEUCHOS_TEST_FOR_EXCEPTION(basisOrder != 1, std::runtime_error, "The order of the basis functions for alphaROM  is not Lagrange(1)");
  // TEUCHOS_TEST_FOR_EXCEPTION(numOfElements != Ru_, std::runtime_error, "The value for Ru for alphaROM does not match that for b");

  // // Create the time DiscreteSpace
  // Array<BasisFamily> time_basisArray(Ru_);
  // for(int i = 0; i < Ru_; i++)
  //   time_basisArray[i] = time_basis_;
  // DiscreteSpace time_DS(timeMesh_, time_basisArray, epetraVecType_);

  // // Establish the values of alphaROM
  // Expr alphaROM = new DiscreteFunction(time_DS, 0.0, alias);
  // // The get() links to the DiscreteFunction's actual vector (shallow)
  // Vector<double> alphaROM_Vec = getDiscreteFunctionVector(alphaROM);
  // for(int i = 0; i < alphaROM_Vec.space().numLocalElements(); i++)
  //   alphaROM_reader >> alphaROM_Vec[i];

  // // Compare alphaROM and alphaOPT
  // double alphaError = L2Norm(timeMesh_, interior_, alphaROM-alphaOPT, quad_);

  // return alphaError;







  
}