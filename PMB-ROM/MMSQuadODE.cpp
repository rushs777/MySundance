#include "MMSQuadODE.hpp"

MMSQuadODE::MMSQuadODE(Teuchos::Array<Expr> phi, Expr uB, Expr q, Expr t, Mesh mesh, bool MatrixAndTensorInFile, int verbosity, int quadOrder) 
    : QuadraticODERHSBase(phi.size(), verbosity),
      interior_(new MaximalCellFilter()),
      phi_(phi),
      uB_(uB),
      q_(q),
      t_(t),
      mesh_(mesh),
      forceIP_(phi.size()),
      MatrixAndTensorInFile_(MatrixAndTensorInFile),
      quad_(new GaussianQuadrature(quadOrder))
  {
    //t_ = new Sundance::Parameter(0.0);
    t_.setParameterValue(0.0);
    Expr x = new CoordExpr(0,"x");
    Expr y = new CoordExpr(1,"y");
    /*double nu = 1.0;
    q_ = List((-36*Cos(y)*Sin(t_)*Sin(x) + (((44 + 30*Cos(t_) + 8*Cos(4*t_) + 30*Cos(5*t_))*
					     Cos(x) + 9*((3 + 2*Cos(t_) + 2*Cos(5*t_))*Cos(3*x) + Cos(5*x)))*
					    Power(Sin(x),3) + 9*Cos(6*t_)*Cos(2*x)*Power(Sin(2*x),3))*Power(Sin(y),2) - 
	       3*(8*nu*Cos(2*t_)*(-1 + 2*Cos(2*x)) + 6*nu*Cos(3*t_)*(-1 + 5*Cos(4*x)) + 
		  8*Sin(2*t_)*Power(Sin(x),2) + 9*Sin(3*t_)*Power(Sin(2*x),2))*Sin(2*y))/36.,
	      (6*nu*(2*Cos(2*t_)*(-1 + 2*Cos(2*y))*Sin(2*x) + 
		     3*Cos(3*t_)*(-4 + 5*Cos(2*y))*Sin(4*x)) - 18*Cos(x)*Sin(t_)*Sin(y) + 
	       3*(4*Sin(2*t_)*Sin(2*x) + 9*Sin(3*t_)*Sin(4*x))*Power(Sin(y),2) + 
	       2*Cos(y)*(Cos(2*t_)*(4*Cos(2*t_) - 3*Cos(3*t_)*(-3 - 6*Cos(2*x) + Cos(4*x)))*
			 Power(Sin(x),2) + 9*Power(Cos(3*t_),2)*Power(Sin(2*x),2))*Power(Sin(y),3))/
	      18.);
    */
    
    // mesh.spatialDim() returns n for nD
    // Define grad operator
    Expr grad = gradient(mesh_.spatialDim());

    for(int r = 0; r < phi_.size(); r++)
      {
	forceIP_[r] = FunctionalEvaluator(mesh_, Integral(interior_, q_*phi_[r], quad_));
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
	    A->setElement(i,r,gradIP(phi_[i], phi_[r]));

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
   * gradIP peforms the IP -(grad*f, grad*g)
   ********************************************************************************/
double MMSQuadODE::gradIP(Expr f, Expr g)
  {
    // mesh.spatialDim() returns n for nD
    int dim = mesh_.spatialDim();

    // Define our differential operators; note Derivative(x=0)
    Expr grad = gradient(dim);

    FunctionalEvaluator IP = FunctionalEvaluator(mesh_, Integral(interior_, -colonProduct(outerProduct(grad,f),outerProduct(grad,g)), quad_));
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
