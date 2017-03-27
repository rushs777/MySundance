#include "MMSQuadODE.hpp"

MMSQuadODE::MMSQuadODE(Teuchos::Array<Expr> phi, Expr uB, Expr q, Expr t, double deltat, Mesh mesh, bool MatrixAndTensorInFile, int verbosity, int quadOrder) 
    : QuadraticODERHSBase(phi.size(), verbosity),
      interior_(new MaximalCellFilter()),
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
    //t_ = new Sundance::Parameter(0.0);
    // NEED TO ADD A NU VALUE
    t_.setParameterValue(0.0);
    tNext_ = new Sundance::Parameter(deltat_);
    Expr x = new CoordExpr(0,"x");
    Expr y = new CoordExpr(1,"y");
    
    // mesh.spatialDim() returns n for nD
    // Define grad operator
    Expr grad = gradient(mesh_.spatialDim());
    
    for(int r = 0; r < phi_.size(); r++)
      {
	Expr integrand = -outerProduct(grad,uB_)*phi_[r]*uB_ + q_*phi_[r] - (1.0)*colonProduct(outerProduct(grad,uB_),outerProduct(grad,phi_[r]));
	forceIP_[r] = FunctionalEvaluator(mesh_, Integral(interior_, integrand, quad_));
      }

  }

Vector<double> MMSQuadODE::evalForceTerm(const double& t) const
  {
    SUNDANCE_ROOT_MSG3(getVerbosity(), "start eval force");
    t_.setParameterValue(t);

    Vector<double> rtn = space().createMember();
    SUNDANCE_ROOT_MSG3(getVerbosity(), "vec size: " << 8*space().dim());
    for(int r = 0; r < phi_.size(); r++)
      {
	rtn[r] = forceIP_[r].evaluate();
      }
    
    SUNDANCE_ROOT_MSG3(getVerbosity(), "end eval force");
    //std::cout << "Here is the value of (vf, phi) " << std::endl << rtn << std::endl;
    return rtn;
  }

void MMSQuadODE::fillMatrixAndTensor(RCP<DenseSerialMatrix>& A, Array<RCP<DenseSerialMatrix> >& T)
  {
    SUNDANCE_ROOT_MSG1(getVerbosity(), "Starting fillMatrixAndTensor");
    // Access the matrix as a DenseSerialMatrix; note that A is a pointer to A_
    // This will create the matrix A = [(grad*phi_i, grad*phi_j)]
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
   ********************************************************************************/
double MMSQuadODE::A_IP(Expr phi_i, Expr phi_j)
  {
    // mesh.spatialDim() returns n for nD
    int dim = mesh_.spatialDim();

    // Define our differential operators; note Derivative(x=0)
    Expr grad = gradient(dim);

    Expr integrand = -outerProduct(grad,phi_j)*phi_i*uB_ - outerProduct(grad,uB_)*phi_i*phi_j - (1.0)*colonProduct(outerProduct(grad,phi_i),outerProduct(grad,phi_j));
    FunctionalEvaluator IP = FunctionalEvaluator(mesh_, Integral(interior_, integrand, quad_));
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
