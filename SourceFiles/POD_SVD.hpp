#ifndef POD_SVD_HPP
#define POD_SVD_HPP

#include "VientoSnapshotIO.hpp" // For writesnap


// Code
#include "PlayaLinearOperatorDecl.hpp"// For LinearOperator
#include "Sundance.hpp"


#include "Teuchos_GlobalMPISession.hpp" //Handles initializing and finalizing MPI. Also defines rcp_dynamic_cast, tuple
#include "PlayaVectorType.hpp" // For VectorType, VectorSpace
#include "PlayaEpetraVectorType.hpp" // For EpetraVectorType
#include "PlayaSerialVectorType.hpp"


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

// Outstanding Questions:
// 1. get_uB() needs to know the number of time steps: one option is to pass it as an argument to
//   the function; the question has actually become where to put this entirely;


/**
 * The class POD_SVD performs a Proper Orthogonal Decomposition via the calculation
 * of the SVD of the matrix B^T S B, where B is a supplied MxR matrix and 
 * S is the RxR mass matrix derived from the supplied DiscreteSpace. 
 * Its purpose is to generate a reduced-order basis.
 */
class POD_SVD
{
public:
  /** Constructor - requires 
   * LinearOperator<double> with an underlying DenseSerialMatrix 
   * DiscreteSpace
   * integer representing the desired level of verbosity
   */
  POD_SVD(const LinearOperator<double> &B, DiscreteSpace &ds, int verbosity = 1);

  /**
   * calculateSVD - performs the SVD of B^T S B using LAPACK's dgesvd function 
   */
  void calculateSVD();

  /**
   * get_lambda - return the singular values of B^T S B; these are also the eigenvalues
   * of B^T S B
   */
  Vector<double> get_lambda() {return lambda_;}

  /**
   * get_singular_vectors() - return the singular vectors as the columns of a matrix
   */
  LinearOperator<double> get_singular_vectors() {return Chi_;}

  /*
   * calculateBasisFunctions - calculate the reduced-order basis functions for the information
   * in B_ by the formula phi_r = (B Chi[r]) / Snorm( B Chi[r]) where Chi[r] denotes the r^{th}
   * column vector of the left singular vector matrix from the SVD
   */
  void calculateBasisFunctions();

  /*
   * get_basis_functions - return the reduced order basis functions produced by
   * calculateBasisFunctions whose relative information content meets the value of tol.
   * Also writes them to fileDir as POD_basis[i]
   */
  Array<Expr> get_basis_functions(double tol, string fileDir);
  


private:
  /** B is the supplied MxR matrix for which the singular values of B^T S B will be calcualted */
  Playa::LinearOperator<double> B_;
  /** The DiscreteSpace for which the mass matrix S should be built */
  DiscreteSpace ds_;
  /** Level of output verbosity; default value is 1 */
  int verbosity_;
  /** S is the RxR mass matrix created by using the supplied 
   *DiscreteSpace to create a LinearProblem*/
  Playa::LinearOperator<double> S_;
  /** phi_ is the Array<Expr> holding the reduced-basis reflecting the information space spanned
   * by B_ */
  Array<Expr> phi_;
  /** lambda_ is the R-length Vector<double> that holds the eigenvalues 
   * (same as the singular values) of B^T S B */
  Vector<double> lambda_;
  /** Chi_ is the MxR matrix that holds the left singular vectors (same as the eigen vectors) */
  Playa::LinearOperator<double> Chi_;

  /** createMassMatrix - builds the mass matrix S_ based off of the supplied DiscreteSpace */
  void createMassMatrix();

  /** Snorm - calculates the weighted mass matrix norm sqrt( v^T S v) by calling Sip.
   * Expects both the matrix S and the Vector v to be based off of Epetra
   */
  double Snorm(const LinearOperator<double> &S, const Vector<double> &v) {return sqrt(Sip(S,v,v));}

  /** Sip - calculates the weighted mass matrix inner product (f,g)_S = f^T S g 
   * Expects both Vectors and the matrix S to be based off of Epetra */
  double Sip(const LinearOperator<double> &S, const Vector<double> &f, const Vector<double> &g);

};

#endif
