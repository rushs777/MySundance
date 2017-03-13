#include "MMSQuadODE.hpp"

MMSQuadODE::MMSQuadODE(Teuchos::Array<Expr> phi, Mesh mesh, bool MatrixAndTensorInFile, int verbosity, int quadOrder) 
    : QuadraticODERHSBase(phi.size(), verbosity), 
      phi_(phi), mesh_(mesh), 
      MatrixAndTensorInFile_(MatrixAndTensorInFile),
      quad_(new GaussianQuadrature(quadOrder))
  {
    t_ = new Sundance::Parameter(0.0);
    Expr x = new CoordExpr(0,"x");
    Expr y = new CoordExpr(1,"y");
    double nu = 1.0;
    q_ = List((-36*cos(y)*sin(t_)*sin(x) + (((44 + 30*cos(t_) + 8*cos(4*t_) + 30*cos(5*t_))*
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
  }

Vector<double> MMSQuadODE::evalForceTerm(const double& t) const
  {
    t_.setParameterValue(t);
    CellFilter interior = new MaximalCellFilter();

    Vector<double> rtn = space().createMember();
    for(int r = 0; r < phi_.size(); r++)
      {
	FunctionalEvaluator IP = FunctionalEvaluator(mesh_, Integral(interior, q_*phi_[r], quad_));		
	rtn[r] = IP.evaluate();
      }

    //std::cout << "Here is the value of (vf, phi) " << std::endl << rtn << std::endl;
    return rtn;
  }

void MMSQuadODE::fillMatrixAndTensor(RCP<DenseSerialMatrix>& A, Array<RCP<DenseSerialMatrix> >& T)
  {
    Out::root() << "Starting fillMatrixAndTensor " << endl;
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
  } // End of fillMatrixAndTensor
	

  /********************************************************************************
   * gradIP peforms the IP (grad*f, grad*g)
   ********************************************************************************/
double MMSQuadODE::gradIP(Expr f, Expr g)
  {

    CellFilter interior = new MaximalCellFilter();
    // mesh.spatialDim() returns n for nD
    int dim = mesh_.spatialDim();

    // Define our differential operators; note Derivative(x=0)
    Expr grad = gradient(dim);

    FunctionalEvaluator IP = FunctionalEvaluator(mesh_, Integral(interior, colonProduct(outerProduct(grad,f),outerProduct(grad,g)), quad_));
    return (IP.evaluate());
  }

  /********************************************************************************
   * tensorIP peforms the IP (f, (h*grad)*g)
   ********************************************************************************/
double MMSQuadODE::tensorIP(Expr f, Expr g, Expr h)
  {
    CellFilter interior = new MaximalCellFilter();
    // mesh.spatialDim() returns n for nD
    // Define grad operator
    Expr grad = gradient(mesh_.spatialDim());

    //      The first one is what KL and I did on 1/27; then I realized the indices were off
    //	FunctionalEvaluator IP = FunctionalEvaluator(mesh, Integral(interior, h*((f*grad)*g),quad));
    FunctionalEvaluator IP = FunctionalEvaluator(mesh_, Integral(interior, f*((h*grad)*g),quad_));
    return (IP.evaluate());
  }
