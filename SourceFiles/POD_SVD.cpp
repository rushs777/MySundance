#include "POD_SVD.hpp" // This classes headerfile





POD_SVD::POD_SVD(const LinearOperator<double> &B, DiscreteSpace &ds, int verbosity)
  : B_(B),
    ds_(ds),
    verbosity_(verbosity)
{
  // Create the mass matrix S_
  createMassMatrix();
}

void POD_SVD::createMassMatrix()
{
  // Make sure that the DiscreteSpace is built off an EpetraVectorType
  const EpetraVectorType* check = dynamic_cast<const EpetraVectorType*>(ds_.vecType().ptr().get());
  TEUCHOS_TEST_FOR_EXCEPTION(check==0, runtime_error,
			  "The DiscreteSpace given is not built off of EpetraVectorType");
  
  SUNDANCE_ROOT_MSG1(verbosity_, "Creating the mass matrix S.........");
  // Create the mass matrix S from the DiscreteSpace ds_

  // Filter subtype MaximalCellFilter selects all cells having dimension
  // equal to the spatial dimension of the ds_.mesh(). 
  Sundance::CellFilter interior = new Sundance::MaximalCellFilter();

  // Test Functions
  Teuchos::Array<Expr> v;
  for(int i = 0; i < ds_.basis().size(); i++)
    {
      v.push_back(new TestFunction(ds_.basis()[i], "v[" + Teuchos::toString(i) + "]"));
    }	
  Sundance::Expr vlist = new Sundance::ListExpr(v);

  // Unknown Functions
  Teuchos::Array<Expr> u;
  for(int i = 0; i < ds_.basis().size(); i++)
    {
      u.push_back(new UnknownFunction(ds_.basis()[i], "u[" + Teuchos::toString(i) + "]"));
    }	
  Sundance::Expr ulist = new Sundance::ListExpr(u);
      
  // Evaluation scheme. The parameter says for what degree of polynomials it will be exact for
  Sundance::QuadratureFamily quad = new Sundance::GaussianQuadrature(4);
  // Define what you want to integrate
  Sundance::Expr integrand = vlist*ulist;
  Sundance::Expr eqn = Integral(interior,integrand,quad);
  // Define Empty BC since I want the mass matrix
  Sundance::Expr bc;
  // Define the problem
  Sundance::LinearProblem prob(ds_.mesh(),eqn,bc,vlist,ulist,ds_.vecType());

  // Get the mass matrix
  S_ = prob.getOperator();

  SUNDANCE_ROOT_MSG2(verbosity_, "The mass matrix S is " << S_.range().dim() << " by " << S_.domain().dim() );
  SUNDANCE_ROOT_MSG2(verbosity_, "Number of vertices in S: " + Teuchos::toString(ds_.mesh().numCells(0)));
  SUNDANCE_ROOT_MSG2(verbosity_, "Number of edges in S: " + Teuchos::toString(ds_.mesh().numCells(1)));
  //SUNDANCE_ROOT_MSG3(verbosity_, S_);
  if(verbosity_>=3)
    cout << S_ << endl;
  SUNDANCE_ROOT_MSG2(verbosity_, "B is " << B_.range().dim() << " by " << B_.domain().dim() );
  //  SUNDANCE_ROOT_MSG3(verbosity_, B_);
  if(verbosity_>=3)
    cout << B_ << endl;
}



void POD_SVD::calculateSVD()
{
  /*
   * The premise here is that instead of solving R eigenvalue problems, we solve one SVD
   * Since B^T S B is symmetric, positive semi-definite, it has all non-negative ews
   * Singular values of B^T S B are the square roots of the ews of 
   * A^T A = (V^T S V)^T (V^T S V) = V^T S^2 V
   * which has eigenvalues that are just the square of the ews of B^T S B
   * Therefore the eigenvalues and singular values of B^T S B are the same
   */

  // Let A = B^T S B; then A = Chi Sigma ChiT
  /* Note: Due to an issue in LAPACK dgesvd, the left and right matrices aren't
   * exactly maintaining the rule that U.V^T = I when only a few of the singular
   * values are noticably greater than 0. Thus we will only be working with
   * U = Chi
   */
  Playa::LinearOperator<double> ChiT;
  SUNDANCE_ROOT_MSG2(verbosity_, "Calculating A = B^T*S*B.........");
  Playa::LinearOperator<double> A = denseMatrixMatrixProduct(B_.transpose(), epetraDenseProduct(S_,B_) ); //denseMatrixMatrixProduct checks the dimensions
  SUNDANCE_ROOT_MSG1(verbosity_, "Getting the SVD of B^T*S*B"); 
  denseSVD(A, Chi_, lambda_, ChiT);

  // Check to make sure that the left singular vectors (same as the eigenvectors) are orthogonal
  LinearOperator<double> checkMtx = denseMatrixMatrixProduct(Chi_.transpose(),Chi_);
  RCP<DenseSerialMatrix> checkMtxPtr = DenseSerialMatrix::getConcretePtr(checkMtx);
  double* allEntries = checkMtxPtr->dataPtr();
  double entry;

  bool orthogonal = true;
  // Remeber that LinearOperator stores the information column wise in an Array
  for(int j = 0; j < checkMtxPtr->numCols(); j++)
    {
      for(int i = 0; i < checkMtxPtr->numRows(); i++)
	{
	  entry = allEntries[i+checkMtxPtr->numRows()*j];
	  if( fabs(entry) > 1e-10 && i != j)
	    {
	      orthogonal = false;
	      Out::os() << "Chi * Chi.transpose() != I " << endl
			<< "Chi[" << i << "]*Chi[" << j << "] = "
			<< entry << " and not 0 " << endl;
	    }
	  else if( fabs(entry - 1.0) > 1e-10 && i == j)
	    {
	      orthogonal = false;
	      Out::os() << "Chi * Chi.transpose() != I " << endl
			<< "Chi[" << i << "]*Chi[" << i << "] = "
			<< entry << " and not 1 " << endl;
	    }
	}
    }

  TEUCHOS_TEST_FOR_EXCEPT(orthogonal==false);
  
  if(verbosity_ >= 3)
    cout << "Here is ChiT.Chi " << endl << checkMtx << endl;
}

void POD_SVD::calculateBasisFunctions()
{
  // Needed for checking that reduced-basis functions are orthogonal wrt L2Norm
  CellFilter interior = new MaximalCellFilter();
  QuadratureFamily quad = new GaussianQuadrature(6);
  
  
  // Find phi_r for r = 1:R (dimension of the mass matrix, rows of B_)
  // Start by getting chi_r from Chi
  Playa::Vector<double> chi_r = Chi_.range().createMember();  //B_.domain().createMember();
  Playa::Vector<double> ej = Chi_.domain().createMember();
  Playa::Vector<double> phiVec = B_.range().createMember();

  Playa::VectorType<double> vecTypeSerial = new Playa::SerialVectorType();
  Playa::VectorSpace<double> dense_S_domain = vecTypeSerial.createEvenlyPartitionedSpace(Playa::MPIComm::world(), S_.domain().dim() );
  phi_.resize(Chi_.domain().dim());

  for(int r = 0; r < Chi_.domain().dim(); r++)
    {
      ej.zero();
      ej[r] = 1.0;

      // Get the r singular vector from Chi_
      chi_r.zero();
      Chi_.apply(ej,chi_r);

      // Now use (10) from "A numerical investigation of velocity-pressure ROM ..."
      phiVec.zero();
      B_.apply(chi_r,phiVec);
      phiVec*=(1.0/ Snorm(S_, serialToEpetra(phiVec) ) );

      // Convert the vector to an Expr
      // DiscreteFunction requires EpetraVectorType
      phi_[r] = new DiscreteFunction(ds_, serialToEpetra(phiVec));

      // Check unit length
      TEUCHOS_TEST_FOR_EXCEPTION( fabs(L2Norm(ds_.mesh(), interior, phi_[r], quad) - 1.0) >= 1.0e-6,
				  runtime_error,
				  "||phi["+Teuchos::toString(r)+"]|| = "
				  + Teuchos::toString(L2Norm(ds_.mesh(), interior, phi_[r], quad))
				  + " != 1");      
    }

  // Check that the reduced-basis functions are orthogonal
  for(int i = 0; i < phi_.length(); i++)
    {
      for(int j = 0; j < phi_.length(); j++)
	{
	  if(i != j)
	    {
	      FunctionalEvaluator dotProduct(ds_.mesh(), Integral(interior, phi_[i]*phi_[j], quad));
	      TEUCHOS_TEST_FOR_EXCEPTION( fabs( dotProduct.evaluate() ) > 1e1,
					  runtime_error,
					  "(phi[" + Teuchos::toString(i) + "], phi["
					  + Teuchos::toString(j) + "]) = "
					  + Teuchos::toString(dotProduct.evaluate())
					  + " != 0");
	    }
	}
    }

}


Array<Expr> POD_SVD::get_basis_functions(double tol, string fileDir)
{
  TEUCHOS_TEST_FOR_EXCEPTION(tol < 0 && tol > 1.0,
			     runtime_error,
			     "The value of tol given as " + Teuchos::toString(tol)
			     + " is not in (0,1]");
  
  Array<Expr> rtn;
  rtn.push_back(phi_[0]);
  double totalInformation = lambda_.norm1();
  double sumInformation = lambda_[0];
  double RIC = sumInformation/totalInformation;

  int i = 1;
  while(RIC < tol)
    {
      rtn.push_back(phi_[i]);
      sumInformation += lambda_[i];
      RIC = sumInformation/totalInformation;
      i++;
    }

  // Write the reduced-order basis functions to file
  int fileError = system( ("rm -fr " + fileDir).c_str() );
  TEUCHOS_TEST_FOR_EXCEPTION( fileError == -1, runtime_error,
			      "Failed to delete " + fileDir
			      + "; can also be triggered if this is the inital run" );
  fileError = system( ("mkdir -p " + fileDir).c_str() );
  TEUCHOS_TEST_FOR_EXCEPTION( fileError == -1, runtime_error,
			      "Failed to create " + fileDir); 
  
  string filePrefix = fileDir + "/POD_basis";
  int Ru = i;
  for(int r = 0; r < Ru; r++)
    writeSnap(filePrefix, r, rtn[r]);

  return rtn;

  
  // for(int i = 0; i < lambda_.dim(); i++)
  //   {
  //     sumInformation += lambda_[i];
  //     RIC = sumInformation/totalInformation;
  //     rtn.push_back(phi_[i]);      
  //     if( RIC >= tol )
  // 	break;
  //   }
  
  
}


 
double POD_SVD::Sip(const LinearOperator<double> &S, const Vector<double> &f, const Vector<double> &g)
{
  Vector<double> Sg = S.range().createMember();
  S.apply(g,Sg);
  double rtn = f*Sg;
  return rtn;
}


	






