#include "PlayaSVD.hpp" // This classes headerfile
#include "Teuchos_GlobalMPISession.hpp" //Handles initializing and finalizing MPI. Also defines rcp_dynamic_cast, tuple
#include "PlayaVectorType.hpp" // For VectorType, VectorSpace
#include "PlayaEpetraVectorType.hpp" // For EpetraVectorType
#include "PlayaSerialVectorType.hpp"
//#include "Sundance.hpp" // For MaximalCellFilter; mentioned in svd.hpp, so this program knows about it

// Not sure about these two
#include "PlayaLinearOperatorImpl.hpp" // Needed to create an empty LinearOperator object

#include "PlayaDenseSerialMatrix.hpp" // For DenseSerialMatrix, denseSVD
#include "PlayaDenseSerialMatrixFactory.hpp" // For DenseSerialMatrixFactory
#include "Teuchos_ParameterXMLFileReader.hpp" //Needed for ParameterXMLFileReader, ParameterList
#include "PlayaDenseMatrixMatrixProduct.hpp"
#include "PlayaSerialEpetraAdapter.hpp" // Allows me to convert between Epetra and Serial Vectors
#include <cstdlib>
using std::cout;
using std::endl;

/**Defines the norm ||v||_s = sqrt(v^T*S*v)
 * Everything will be expected to be based of the EpetraVectorType
 */
double Snorm(const LinearOperator<double> &S, const Vector<double> &v)
{
  Vector<double> Sv = S.range().createMember();
  S.apply(v,Sv);
  double rtn = sqrt(v*Sv);
  return rtn;
}

/**Defines the IP (f,g)_s = f^T*S*g
 * Everything will be expected to be based of the EpetraVectorType
 */
double Sip(const LinearOperator<double> &S, const Vector<double> &f, const Vector<double> &g)
{
  Vector<double> Sg = S.range().createMember();
  S.apply(g,Sg);
  double rtn = f*Sg;
  return rtn;
}


namespace Playa
{
  /*
   *
   * S: mass matrix; n x n
   * W: snapshot matrix; n x m
   * lambda: an array to hold the eigenvalues
   * Alpha: a matrix that holds time coefficient vectors of the 
   *        eigenmodes of the POD as columns
   * Phi: a matrix that holds the eigenmodes (eigenvectors) of the POD as columns
   *
   *
   */
  void POD(const LinearOperator<double> &W, Vector<double> &lambda, LinearOperator<double> &Alpha, LinearOperator<double> &Phi, Sundance::Mesh &mesh, int debug)
  {
    SUNDANCE_ROOT_MSG1(debug, "Creating the mass matrix S.........");
    // Create the mass matrix S
    Playa::VectorType<double> vecTypeEpetra = new Playa::EpetraVectorType();
    // Filter subtype MaximalCellFilter selects all cells having dimension equal to the spatial dimension of the mesh. 
    Sundance::CellFilter interior = new Sundance::MaximalCellFilter();
    // BasisFamily used to express our solutions
    //	Sundance::BasisFamily basis = new Sundance::Lagrange(2); // 2nd order PWQL (Piece-Wise Quadratic Lagrange)
    Array<Sundance::BasisFamily> basis; 
    if(mesh.spatialDim()==2)
      basis = List(new Sundance::Lagrange(2), new Sundance::Lagrange(2)); // 2nd order PWQL (Piece-Wise Quadratic Lagrange)
    else if(mesh.spatialDim()==3)
      {
	basis = List(new Sundance::Lagrange(2), new Sundance::Lagrange(2), new Sundance::Lagrange(2)); // 2nd order PWQL (Piece-Wise Quadratic Lagrange)
      }
    else
      {
	//basis = List(new Sundance::Lagrange(1), new Sundance::Lagrange(1));
	basis.push_back(new Sundance::Lagrange(1));
      }

    // Test Functions
    int numTest = mesh.spatialDim();
    Teuchos::Array<Expr> v(numTest);
    for(int i=0; i<v.size(); i++)
      v[i] = new TestFunction(basis[i], "v[" + Teuchos::toString(i) + "]");

	
    Sundance::Expr vlist = new Sundance::ListExpr(v);

    // Unknown Functions
    int numUnknown = mesh.spatialDim();
    Teuchos::Array<Expr> u(numUnknown);
    for(int i=0; i<u.size(); i++)
      u[i] = new UnknownFunction(basis[i], "u[" + Teuchos::toString(i) + "]");

	
    Sundance::Expr ulist = new Sundance::ListExpr(u);
      
    // Evaluation scheme. The parameter says for what degree of polynomials it will be exact for
    Sundance::QuadratureFamily quad = new Sundance::GaussianQuadrature(4);
    // Define what you want to integrate
    Sundance::Expr integrand = vlist*ulist;
    Sundance::Expr eqn = Integral(interior,integrand,quad);
    // Define Empty BC
    Sundance::Expr bc;// since I want the mass matrix
    // Define the problem
    Sundance::LinearProblem prob(mesh,eqn,bc,vlist,ulist,vecTypeEpetra);

    // Get the mass matrix
    Playa::LinearOperator<double> S = prob.getOperator();

    SUNDANCE_ROOT_MSG2(debug, "S is " << S.range().dim() << " by " << S.domain().dim() );
    SUNDANCE_ROOT_MSG2(debug, "W is " << W.range().dim() << " by " << W.domain().dim() );
    std::cout << "Number of vertices in S: " << mesh.numCells(0) << std::endl
	      << "Number of edges in S: " << mesh.numCells(1) << std::endl;

    //	std::cout << "S = " << std::endl << S << std::endl;
    /* Check to make sure I was getting the right mass matrix
       Playa::Vector<double> x = S2.domain().createMember();
       x.randomize();
       Playa::Vector<double> z = S*x - S2*x;
       double errorz = z.norm2();
       std::cout << "Here is ||S*x-S2*x|| = " << errorz << std::endl;
    */
    /*
     * m: the number of time steps; will be the number of columns of W
     * n: the number of grid points; will be the number of rows of W
     */
    //int m = W.domain().dim();
    //int n = W.range().dim();

    //Teuchos::ParameterXMLFileReader xmlReader("anasazi-ml.xml");
    //Teuchos::ParameterList solverParams = xmlReader.getParameters().sublist("Eigensolver");

    Playa::LinearOperator<double> Alphat;
    SUNDANCE_ROOT_MSG2(debug, "Calculating A = W^t*S*W");
    Playa::LinearOperator<double> A = denseMatrixMatrixProduct(W.transpose(), epetraDenseProduct(S,W) ); //denseMatrixMatrixProduct checks the dimensions
    SUNDANCE_ROOT_MSG1(debug, "Getting the dense SVD"); 
    denseSVD(A, Alpha, lambda, Alphat);

    // Find phi_r for r = 1:R
    // Start by getting alpha_r from Alpha
    Playa::Vector<double> alpha_r = W.domain().createMember();
    Playa::Vector<double> ej = Alpha.domain().createMember();
    //Teuchos::RCP<Playa::DenseSerialMatrix> PhiPtr = 				Teuchos::rcp_dynamic_cast<Playa::DenseSerialMatrix>(Phi.ptr());

    Playa::VectorType<double> vecTypeSerial = new Playa::SerialVectorType();
    Playa::VectorSpace<double> dense_S_domain = vecTypeSerial.createEvenlyPartitionedSpace(Playa::MPIComm::world(), S.domain().dim() );


    Playa::DenseSerialMatrixFactory PhiMF(Alpha.domain(), dense_S_domain ); //Alpha.domain().dim() = R
    Phi = PhiMF.createMatrix();
    RCP<DenseSerialMatrix> PhiPtr = DenseSerialMatrix::getConcretePtr(Phi);
    //Teuchos::RCP<Playa::DenseSerialMatrix> PhiPtr = Teuchos::rcp_dynamic_cast<Playa::DenseSerialMatrix>(Phi.ptr());
    //DenseSerialMatrix* PhiPtr = dynamic_cast<DenseSerialMatrix*>(Phi.ptr().get());
    //	Playa::Vector<double> phi_r = dense_S_domain.createMember();
    Playa::Vector<double> phi_r = Phi.range().createMember();
    double* PhiData = PhiPtr->dataPtr();

    for(int count = 0; count < PhiPtr->numCols(); count++)
      {
	ej.zero();
	alpha_r.zero();
	phi_r.zero();
	ej[count]=1.0;

	/*std::cout << "Size of ej " << ej.dim() << std::endl
	  << "Size of alpha_r " << alpha_r.dim() << std::endl
	  << "Size of phi_r " << phi_r.dim() << std::endl;*/

	//std::cout << "Here is e_" << count << std::endl << ej << std::endl;
	Alpha.apply(ej,alpha_r);
	//std::cout << "Here is alpha_" << count << std::endl << alpha_r << std::endl;

	// Now use (10) from "A numerical investigation of velocity-pressure ROM ..."
	W.apply(alpha_r,phi_r);
	//std::cout << "Here is W*alpha_" << count << std::endl << phi_r << std::endl;
	//		phi_r*=(1.0/ ( sqrt(lambda[count])*sqrt(alpha_r*alpha_r) ));
	phi_r*=(1.0/ Snorm(S, serialToEpetra(phi_r) ) );
	//std::cout << "Here is phi_" << count << std::endl << phi_r << std::endl;
	
	// Store phi_r as the r^th column of Phi
	for(int row = 0; row<PhiPtr->numRows(); row++)
	  PhiData[row+PhiPtr->numRows()*count] = phi_r[row];
	
      }	

    // Set a flag to make sure that all of the phi_r are orthonormal
    bool orthonormal = true;
    Vector<double> ei = Phi.domain().createMember();	
    Vector<double> phi_i = Phi.range().createMember();
    Vector<double> phi_j = Phi.range().createMember();
    cout << "Phi is " << PhiPtr->numRows() << "x" << PhiPtr->numCols() << endl;
    cout << "ei is " << ei << endl;
    cout << "ej is " << ej << endl;
    cout << "phi_i is " << phi_i << endl;
    cout << "phi_j is " << phi_j << endl;
    Vector<double> Sx = S.range().createMember();
    for(int i = 0; i<PhiPtr->numCols(); i++)
      {
	ei.zero();
	ei[i] = 1.0;
	phi_i.zero(); // Cause of previous error; apply adds out to that matrix multiply
	Phi.apply(ei,phi_i);

	for(int j = 0; j<PhiPtr->numCols(); j++)
	  {
	    ei.zero();
	    ei[j] = 1.0;
	    phi_j.zero();
	    Phi.apply(ei,phi_j);
	    if( fabs( Sip(S,serialToEpetra(phi_i), serialToEpetra(phi_j)) ) < 1.0e-10 )
	      {
		cout << "(phi[" << i << "], phi[" << j << "])_S = " << Sip(S,serialToEpetra(phi_i), serialToEpetra(phi_j)) << endl;
		orthonormal = false;
	      }
	  }

	if( fabs(Snorm(S, serialToEpetra(phi_i) ) -1.0) > 1.0e-10)
	  {
	    cout << "phi[" << i << "].norm2() = " << phi_i.norm2() << endl;
	    orthonormal = false;
	  }

	if(i>-1)
	  {
	    Vector<double> Sphi_r = S.range().createMember();
	    S.apply(serialToEpetra(phi_i), Sphi_r);
	    Vector<double> Wt_Sphi_r = W.domain().createMember();
	    Vector<double> result = W.range().createMember();
	    W.applyTranspose(epetraToSerial(Sphi_r), Wt_Sphi_r);
	    W.apply(Wt_Sphi_r, result);

	    if( (result-lambda[i]*phi_i).norm2() > 1.0e-4)
	      {
		cout << "Are my coefficient vectors solving equation 8? " << endl
		     << "Testing for r = " << i << endl;
			
		cout << "Coefficient vector for r = " << i << " does not satisfy equation 8 "
		     << endl << "||lhs - rhs||_2 = " << (result-lambda[i]*phi_i).norm2() << endl;
	      }
	  }
      }

    std::cout << "Is everything orthnormal? " << orthonormal << std::endl;
    /*std::cout << "Phi is " << PhiPtr->numRows() << "x" << PhiPtr->numCols() << std::endl;
    std::cout << "Alpha is " << Alpha.range().dim() << "x" << Alpha.domain().dim() << std::endl;
    std::cout << "W is " << W.range().dim() << "x" << W.domain().dim() << std::endl;
    std::cout << "lambda is " << lambda.dim() << "x1" << std::endl;*/






  } //End of POD
	


} // End of File



	






