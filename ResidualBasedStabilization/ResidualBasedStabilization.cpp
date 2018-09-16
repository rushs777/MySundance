#include "ResidualBasedStabilization.hpp"

ResidualBasedStabilization::ResidualBasedStabilization(velocityROM velROM, PMB pressROM, Expr forceTerm, Expr t, int K, int verbosity)
  : velROM_(velROM),
    pressROM_(pressROM),
    forceTerm_(forceTerm),
    t_(t),
    deltat_(deltat),
    K_(K),
    verbosity_(verbosity)
{
  t_.setParameterValue(0.0);
}

void ResidualBasedStabilization::initialize()
{
  Array<Expr> uRO_(velROM_.get_uRO() );
  Array<Expr> pRO_(pressROM_.get_pRO() );
  DiscreteSpace RmDS_ = velROM_.get_ds();
  DiscreteSpace RcDS_ = pressROM_.get_ds();
  double nu_ = 1.0; //Need to fix this

  // Create the snapshot matrix for the momentum residual
  Array<Vector<double> > alpha(velROM_.get_alpha());
  int nSteps = alpha.length();
  int numNodes = RmDS_.mesh().numCells(0) + RMDS_.mesh().numCells(1);
  int spatialDim = RmDS_.mesh().spatialDim();
  VectorType<double> vecTypeSerial = new SerialVectorType();
  VectorSpace<double> domain = vecTypeSerial.createEvenlyPartitionedSpace(MPIComm::world(), nSteps);

  VectorSpace<double> Rm_range = vecTypeSerial.createEvenlyPartitionedSpace(MPIComm::world(), numNodes*spatialDim);
  RCP<Playa::MatrixFactory<double> > Rm_mf = vecTypeSerial.createMatrixFactory(domain,Rm_range);
  LinearOperator<double> momentumSnapshotMtx = Rm_mf->createMatrix();
  RCP<DenseSerialMatrix> momentumPtr = DenseSerialMatrix::getConcretePtr(momentumSnapshotMtx);

  for(int time = 0; time < nSteps; time++)
    {
      Vector<double> column = calculateRm(time);

      for(int i = 0; i < momentumSnapshotMtx.range().dim(); i++)
	{
	  momentumPtr->setElement(i,time,column[i]);
	}
    }

  // Build the POD modes for the residual
  Playa::LinearOperator<double> Rm_U;
  Playa::LinearOperator<double> Zeta;
  Playa::Vector<double> Rm_sigma;

  // Data matrix and DiscreteSpace need to be defined
  SUNDANCE_ROOT_MSG1(verbosity_, "Entering POD for momentum residual");
  POD(momentumSnapshotMtx,Rm_sigma,Rm_U,Zeta,RmDS_,verbosity_);
  SUNDANCE_ROOT_MSG2(verbosity_, "POD finished for momentum residual");

  Vector<double> ej = Zeta.domain().createMember();
  Vector<double> zetaCoeff = Zeta.range().createMember();
  zeta_.resize(K); // These are the POD basis functions


  // Get the Expr zeta_r(x)
  CellFilter interior = new MaximalCellFilter();
  QuadratureFamily quad = new GaussianQuadrature(6);
  for(int r = 0; r<K_; r++)
    {
      ej.zero();
      ej[r] = 1.0;
      zetaCoeff.zero();
      Zeta.apply(ej,zetaCoeff);
      zeta_[r] = new DiscreteFunction(RmDS_, serialToEpetra(zetaCoeff)); //DiscreteFunction requires Epetra vectors
      TEUCHOS_TEST_FOR_EXCEPTION( fabs(L2Norm(RmDS_.mesh(), interior, zeta_[r], quad) - 1.0) >= 1.0e-6,
				  runtime_error, "||zeta["+Teuchos::toString(r)+"]|| = " + Teuchos::toString(L2Norm(RmDS_.mesh(), interior, zeta_[r], quad)) + " != 1");
    }



  // Create the snapshot matrix for the continuity residual
  VectorSpace<double> Rc_range = vecTypeSerial.createEvenlyPartitionedSpace(MPIComm::world(), numNodes);
  RCP<Playa::MatrixFactory<double> > Rc_mf = vecTypeSerial.createMatrixFactory(domain,Rc_range);
  LinearOperator<double> continuitySnapshotMtx = Rc_mf->createMatrix();
  RCP<DenseSerialMatrix> continuityPtr = DenseSerialMatrix::getConcretePtr(continuitySnapshotMtx);

  for(int time = 0; time < nSteps; time++)
    {
      Vector<double> column = calculateRc(time);

      for(int i = 0; i < continuitySnapshotMtx.range().dim(); i++)
	{
	  continuityPtr->setElement(i,time,column[i]);
	}
    }


  // Build the POD modes for the residual
  Playa::LinearOperator<double> Rc_U;
  Playa::LinearOperator<double> Eta;
  Playa::Vector<double> Rc_sigma;

  // Data matrix and DiscreteSpace need to be defined
  SUNDANCE_ROOT_MSG1(verbosity_, "Entering POD for continuity residual");
  POD(continuitySnapshotMtx,Rc_sigma,Rc_U,Eta,RcDS_,verbosity_);
  SUNDANCE_ROOT_MSG2(verbosity_, "POD finished for continuity residual");

  Vector<double> ek = Eta.domain().createMember();
  Vector<double> etaCoeff = Eta.range().createMember();
  eta_.resize(K); // These are the POD basis functions

  // Get the Expr eta_r(x)
  for(int r = 0; r<K_; r++)
    {
      ek.zero();
      ek[r] = 1.0;
      etaCoeff.zero();
      Eta.apply(ek,etaCoeff);
      eta_[r] = new DiscreteFunction(RcDS_, serialToEpetra(zetaCoeff)); //DiscreteFunction requires Epetra vectors
      TEUCHOS_TEST_FOR_EXCEPTION( fabs(L2Norm(RcDS_.mesh(), interior, zeta_[r], quad) - 1.0) >= 1.0e-6,
				  runtime_error, "||zeta["+Teuchos::toString(r)+"]|| = " + Teuchos::toString(L2Norm(RcDS_.mesh(), interior, zeta_[r], quad)) + " != 1");
    }

}
    

void ResidualBasedStabilization::generateBasisFunctions()
{
  velBasisFuncs_ = velROM_.get_phi();

  for(int i = 0; i < K_; i++)
    velBasisFuncs.push_back(zeta_[i]);

  pressBasisFuncs_ = pressROM_.get_psi();

  for(int i = 0; i< K_; i++)
    pressBasisFuncs.push_back(eta_[i]);

  velBasisFuncs_ = GramSchmidt(velBasisFuncs_, RmDS_.mesh(), 6);
  pressBasisFuncs_ = GramSchmidt(pressBasisFuncs_, RcDS_.mesh(),6);
}



Array<Expr> get_uRO()
{
  Array<Expr> stabilized_uRO = velROM_.get_uRO();
  Array<Vector<double> > alpha = velROM_.get_alpha();
  int Nr = velROM_.get_phi().length() - K_;
  
  for(int time = 0; time < stabilized_uRO.length(); time++)
    {
      for(int r=Nr; r<Nr+K_; r++)
	{
	  stabilized_uRO[time] = stabilized_uRO[time] + alpha[time][r]*????
	}
    }
}


Vector<double> ResidualBasedStabilization::calculateRm(int timeIndex)
{
  Expr grad = gradient(RmDS_.mesh().spatialDim());
  
  Expr rtnExpr = (uRO_[timeIndex+1] - uRO_[timeIndex])/deltat_ + outerProduct(grad, uRO_[timeIndex])*uRO_[timeIndex] + (grad*pRO_[timeIndex]) - nu_*( (grad*div(uRO_[timeIndex])) - curl(curl(uRO_[timeIndex])) );

  L2Projector projector(RmDS_, rtnExpr);
  Expr rtnProj = projector.project();
  Vector<double> rtnVec = getDiscreteFunctionVector(rtnProj);    
}

Vector<double> ResidualBasedStabilization::calculateRc(int timeIndex)
{
  Expr rtnExpr = div(uRO_[timeIndex]);

  L2Projector projector(RcDS_, rtnExpr);
  Expr rtnProj = projector.project();
  Vector<double> rtnVec = getDiscreteFunctionVector(rtnProj);
}


Array<Expr> ResidualBasedStabilization::GramSchmidt(Array<Expr> v, Mesh mesh, int quadOrder)
{
  if(v.length()<1)
    {
      Out::root() << "Array of length 0 passed to GramSchmidt" << endl;
      return v;
    }

  CellFilter interior = new MaximalCellFilter();
  QuadratureFamily quad = new GaussianQuadrature(quadOrder);

  v[0] = v[0]*(1.0/L2Norm(mesh, interior, v[0], quad));

  for(int i=1; i<v.length(); i++)
    {
      // make v[i] orthogonal
      for(int j=0; j<i; j++)
	{
	  FunctionalEvaluator IP = FunctionalEvaluator(mesh, Integral(interior, v[i]*v[j], quad));
	  v[i] = v[i] - IP.evaluate()*v[j];
	}

      // make v[i] orthonormal
      v[i] = v[i]*(1.0/L2Norm(mesh, interior, v[i], quad));
    }

  return v;
}
