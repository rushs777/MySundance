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
  double nu_ = 1.0; //Need to fix this

  // Create the snapshot matrix for the momentum residual
  Array<Vector<double> > alpha(velROM_.get_alpha());
  int nSteps = alpha.length();
  int numNodes = velROM_.get_ds().mesh().numCells(0) + velROM_.get_ds().mesh().numCells(1);
  int spatialDim = velROM_.get_ds().mesh().spatialDim();
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
  POD(momentumSnapshotMtx,Rm_sigma,Rm_U,Zeta,velROM_.get_ds(),verbosity_);
  SUNDANCE_ROOT_MSG2(verbosity_, "POD finished for momentum residual");

  Vector<double> ej = Zeta.domain().createMember();
  Vector<double> zetaCoeff = Zeta.range().createMember();
  zeta_.resize(K); // These are the POD basis functions
  
  
  
  
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
  POD(continuitySnapshotMtx,Rc_sigma,Rc_U,Eta,pressROM_.get_ds(),verbosity_);
  SUNDANCE_ROOT_MSG2(verbosity_, "POD finished for continuity residual");

  Vector<double> ek = Eta.domain().createMember();
  Vector<double> etaCoeff = Eta.range().createMember();
  eta_.resize(K); // These are the POD basis functions

}
    









Vector<double> ResidualBasedStabilization::calculateRm(int timeIndex)
{
  Expr grad = gradient(velROM_.get_ds().mesh().spatialDim());
  
  Expr rtnExpr = (uRO_[timeIndex+1] - uRO_[timeIndex])/deltat_ + outerProduct(grad, uRO_[timeIndex])*uRO_[timeIndex] + (grad*pRO_[timeIndex]) - nu_*( (grad*div(uRO_[timeIndex])) - curl(curl(uRO_[timeIndex])) );

  L2Projector projector(velROM_.get_ds(), rtnExpr);
  Expr rtnProj = projector.project();
  Vector<double> rtnVec = getDiscreteFunctionVector(rtnProj);    
}

Vector<double> ResidualBasedStabilization::calculateRc(int timeIndex)
{
  Expr rtnExpr = div(uRO_[timeIndex]);

  L2Projector projector(pressROM_.get_ds(), rtnExpr);
  Expr rtnProj = projector.project();
  Vector<double> rtnVec = getDiscreteFunctionVector(rtnProj);
}
