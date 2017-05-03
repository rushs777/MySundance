#include "MMSQuadODE.hpp"

/*******************************************************************************************
 *
 * NOTE: IF YOUR CODE IS NOT WORKING MAKE SURE IT IS NOT BECAUSE uB DEPENDS ON TIME
 * CURRENTLY uB = ubar, so the time derivative is zero and not implemented
 *
 ******************************************************************************************/

MMSQuadODE::MMSQuadODE(Teuchos::Array<Expr> phi, Expr uB, Expr q, Expr t, double deltat, Mesh mesh, bool MatrixAndTensorInFile, int verbosity, int quadOrder) 
    : QuadraticODERHSBase(phi.size(), verbosity),
      interior_(new MaximalCellFilter()),
      boundary_(new BoundaryCellFilter()),
      phi_(phi),
      uB_(uB),
      q_(q),
      t_(t),
      deltat_(deltat),
      mesh_(mesh),
      forceIP_(phi.size()),
      MatrixAndTensorInFile_(MatrixAndTensorInFile),
      quad_(new GaussianQuadrature(quadOrder))
  {
    // NEED TO ADD A NU VALUE
    nu_ = 1.0;
    t_.setParameterValue(0.0);
    tNext_ = new Sundance::Parameter(deltat_);
    
    // mesh.spatialDim() returns n for nD
    // Define grad operator
    Expr grad = gradient(mesh_.spatialDim());
    Expr nHat = CellNormalExpr(mesh_.spatialDim(), "nHat");
    
    for(int r = 0; r < phi_.size(); r++)
      {
	//Expr integrand = -outerProduct(grad,uB_)*phi_[r]*uB_ + q_*phi_[r] - nu_*colonProduct(outerProduct(grad,uB_),outerProduct(grad,phi_[r]));
	//	Expr integrand = -(uB_*grad)*uB_*phi_[r] + q_*phi_[r] - nu_*colonProduct(outerProduct(grad,uB_),outerProduct(grad,phi_[r]));
	Expr integrand_interior = -(uB_*grad)*uB_*phi_[r] + q_*phi_[r] - nu_*colonProduct(outerProduct(grad,uB_),outerProduct(grad,phi_[r]));
	Expr integrand_boundary = nu_*(nHat*((phi_[r]*grad)*uB_));
	//integrand_boundary = 0.0;

	
	//	forceIP_[r] = FunctionalEvaluator(mesh_, Integral(interior_, integrand, quad_));
	forceIP_[r] = FunctionalEvaluator(mesh_, Integral(interior_, integrand_interior, quad_) + Integral(boundary_, integrand_boundary, quad_));
      }

  }

Vector<double> MMSQuadODE::evalForceTerm(const double& t) const
  {
    //SUNDANCE_ROOT_MSG3(getVerbosity(), "start eval force");
    t_.setParameterValue(t);

    Vector<double> rtn = space().createMember();
    //SUNDANCE_ROOT_MSG3(getVerbosity(), "vec size: " << 8*space().dim());
    for(int r = 0; r < phi_.size(); r++)
      {
	rtn[r] = forceIP_[r].evaluate();
      }
    
    //SUNDANCE_ROOT_MSG3(getVerbosity(), "end eval force");
    //std::cout << "Here is the value of (vf, phi) " << std::endl << rtn << std::endl;
    return rtn;
  }

void MMSQuadODE::fillMatrixAndTensor(RCP<DenseSerialMatrix>& A, Array<RCP<DenseSerialMatrix> >& T)
  {
    SUNDANCE_ROOT_MSG1(getVerbosity(), "Starting fillMatrixAndTensor");
    // Access the matrix as a DenseSerialMatrix; note that A is a pointer to A_
    string fileDir = "A_and_T";
    system( ("mkdir -p " + fileDir).c_str() );

    if(!MatrixAndTensorInFile_)
      {
	SUNDANCE_ROOT_MSG1(getVerbosity(), "Creating A");
	for(int i = 0; i<A->numRows(); i++)
	  for(int r = 0; r<A->numCols(); r++)
	    A->setElement(i,r,A_IP(phi_[i], phi_[r])); //HOW TO HANDLE PASSING TIME?

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
	

  /********************************************************************************
   * A_IP peforms the IP -(phi_i, uB*(grad*phi_j)) - (phi_i, phi_j*(grad*uB_)) - nu*(grad*f, grad*g)
   * Does not have the time derivative component implemented
   ********************************************************************************/
double MMSQuadODE::A_IP(Expr phi_i, Expr phi_j)
  {
    // mesh.spatialDim() returns n for nD
    int dim = mesh_.spatialDim();

    // Define our differential operators; note Derivative(x=0)
    Expr grad = gradient(dim);
    Expr nHat = CellNormalExpr(dim, "nHat");

    //Expr integrand = -outerProduct(grad,phi_j)*phi_i*uB_ - outerProduct(grad,uB_)*phi_i*phi_j - nu_*colonProduct(outerProduct(grad,phi_i),outerProduct(grad,phi_j));
    Expr integrand_interior = -phi_i*((uB_*grad)*phi_j) - phi_i*((phi_j*grad)*uB_) - nu_*colonProduct(outerProduct(grad,phi_i),outerProduct(grad,phi_j));
    Expr integrand_boundary = nu_*(nHat*((phi_i*grad)*phi_j));
    //integrand_boundary = 0.0;
    
    FunctionalEvaluator IP = FunctionalEvaluator(mesh_, Integral(interior_, integrand_interior, quad_) + Integral(boundary_, integrand_boundary, quad_));
    return (IP.evaluate());
  }

  /********************************************************************************
   * tensorIP peforms the IP -(f, (h*grad)*g)
   ********************************************************************************/
double MMSQuadODE::tensorIP(Expr f, Expr g, Expr h)
  {
    // mesh.spatialDim() returns n for nD
    // Define grad operator
    Expr grad = gradient(mesh_.spatialDim());

    FunctionalEvaluator IP = FunctionalEvaluator(mesh_, Integral(interior_, -f*((h*grad)*g),quad_));
    return (IP.evaluate());
  }
