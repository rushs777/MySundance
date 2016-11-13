#include "PlayaSVD.hpp" // This classes headerfile
#include "Teuchos_GlobalMPISession.hpp" //Handles initializing and finalizing MPI. Also defines rcp_dynamic_cast, tuple
#include "PlayaVectorType.hpp" // For VectorType, VectorSpace
#include "PlayaSerialVectorType.hpp" // For SerialVectorType
//#include "Sundance.hpp" // For MaximalCellFilter; mentioned in svd.hpp, so this program knows about it

// Not sure about these two
#include "PlayaLinearOperatorImpl.hpp" // Needed to create an empty LinearOperator object

#include "PlayaDenseSerialMatrix.hpp" // For DenseSerialMatrix, denseSVD
#include "PlayaDenseSerialMatrixFactory.hpp" // For DenseSerialMatrixFactory
#include "Teuchos_ParameterXMLFileReader.hpp" //Needed for ParameterXMLFileReader, ParameterList
#include "PlayaDenseMatrixMatrixProduct.hpp"
//#include "/home/sirush/PhDResearch/DenseMatrixMatrixProduct/PlayaDenseMatrixMatrixProduct.hpp"
#include <cstdlib>

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
void POD(const LinearOperator<double> &W, Vector<double> &lambda, LinearOperator<double> &Alpha, LinearOperator<double> &Phi, Sundance::Mesh &mesh)
{
	// Create the mass matrix S
	Playa::VectorType<double> vecType = new Playa::SerialVectorType();
	// Filter subtype MaximalCellFilter selects all cells having dimension equal to the spatial dimension of the mesh. 
      	Sundance::CellFilter interior = new Sundance::MaximalCellFilter();
      	// BasisFamily used to express our solutions
      	Sundance::BasisFamily basis = new Sundance::Lagrange(2); // 2nd order PWQL (Piece-Wise Quadratic Lagrange)

	// Test Functions
	int numTest = 1;
	Teuchos::Array<Expr> v(numTest);
	for(int i=0; i<v.size(); i++)
		v[i] = new TestFunction(basis, "v[" + Teuchos::toString(i) + "]");

	
	Sundance::Expr vlist = new Sundance::ListExpr(v);

	// Unknown Functions
	int numUnknown = 1;
	Teuchos::Array<Expr> u(numUnknown);
	for(int i=0; i<u.size(); i++)
		u[i] = new UnknownFunction(basis, "u[" + Teuchos::toString(i) + "]");

	
	Sundance::Expr ulist = new Sundance::ListExpr(u);
      
	// Evaluation scheme. The parameter says for what degree of polynomials it will be exact for
	Sundance::QuadratureFamily quad = new Sundance::GaussianQuadrature(4);
	// Define what you want to integrate
	Sundance::Expr integrand = vlist*ulist;
	Sundance::Expr eqn = Integral(interior,integrand,quad);
	// Define Empty BC
	Sundance::Expr bc;// since I want the mass matrix
	// Define the problem
	Sundance::LinearProblem prob(mesh,eqn,bc,vlist,ulist,vecType);
	// Get the mass matrix
	Playa::LinearOperator<double> S = prob.getOperator();
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

	std::cout << "Hi from POD() " << std::endl;

	Playa::LinearOperator<double> Alphat;
	Playa::LinearOperator<double> A = denseMatrixMatrixProduct( denseMatrixMatrixProduct(W.transpose(),S), W); //denseMatrixMatrixProduct checks the dimensions
	denseSVD(A, Alpha, lambda, Alphat);

	// Find phi_r for r = 1:R
	// Start by getting alpha_r from Alpha
	Playa::Vector<double> alpha_r = W.domain().createMember();
	Playa::Vector<double> ej = Alpha.domain().createMember();
	//Teuchos::RCP<Playa::DenseSerialMatrix> PhiPtr = 				Teuchos::rcp_dynamic_cast<Playa::DenseSerialMatrix>(Phi.ptr());


	Playa::DenseSerialMatrixFactory PhiMF(Alpha.domain(), S.domain() ); //Alpha.domain().dim() = R
	Phi = PhiMF.createMatrix();
	DenseSerialMatrix* PhiPtr = dynamic_cast<DenseSerialMatrix*>(Phi.ptr().get());
	Playa::Vector<double> phi_r = S.domain().createMember();
	double* PhiData = PhiPtr->dataPtr();

	for(int count = 0; count < ej.dim(); count++)
{
	ej.zero();
	alpha_r.zero();
	phi_r.zero();
	ej[count]=1.0;
	//std::cout << "Here is e_" << count << std::endl << ej << std::endl;
	Alpha.apply(ej,alpha_r);
	//std::cout << "Here is alpha_" << count << std::endl << alpha_r << std::endl;

	// Now use (10) from "A numerical investigation of velocity-pressure ROM ..."
	W.apply(alpha_r,phi_r);
	//std::cout << "Here is W*alpha_" << count << std::endl << phi_r << std::endl;
	phi_r *= (1.0/phi_r.norm2());
	//std::cout << "Here is phi_" << count << std::endl << phi_r << std::endl;

	// Store phi_r as the r^th column of Phi
	for(int row = 0; row<PhiPtr->numRows(); row++)
		PhiData[row+PhiPtr->numRows()*count] = phi_r[row];

}

	/*std::cout << "Phi is " << PhiPtr->numRows() << "x" << PhiPtr->numCols() << std::endl;
	std::cout << "Alpha is " << Alpha.range().dim() << "x" << Alpha.domain().dim() << std::endl;
	std::cout << "W is " << W.range().dim() << "x" << W.domain().dim() << std::endl;
	std::cout << "lambda is " << lambda.dim() << "x1" << std::endl;*/
	//PhiData[9+10*8] = 1.0;







} //End of POD
	


} // End of File



	






