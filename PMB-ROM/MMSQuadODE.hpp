#ifndef MMSQuadODE_2D_HPP
#define MMSQuadODE_2D_HPP

#include "PlayaSerialVectorType.hpp"
#include "PlayaDenseSerialMatrix.hpp"
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaSerialEpetraAdapter.hpp" // Allows me to convert between Epetra and Serial Vectors
#include "PlayaOut.hpp"
#include <vector>


#include "denseSerialMatrixIO.hpp"
#include "QuadraticODERHSBase.hpp"
#include "MathematicaConverter.hpp"

/** */
class MMSQuadODE : public QuadraticODERHSBase
{
public:
  /**Constructor*/
  MMSQuadODE(Teuchos::Array<Expr> phi, Expr uB, Expr q, Expr t, double deltat, Mesh mesh, bool MatrixAndTensorInFile = false, int verbosity = 1, int quadOrder = 6);

  /**Project the force term onto the basis*/
  Vector<double> evalForceTerm(const double& t) const;

  /**Create A and T according to this problem*/
  void fillMatrixAndTensor(RCP<DenseSerialMatrix>& A, Array<RCP<DenseSerialMatrix> >& T);

private:
  /********************************************************************************
   *
   * Teuchos::Array<Expr> phi - Array of basis function obtained from a POD;
   *                            Each Expr has an underlying DiscreteFunction
   * Mesh mesh - Considered the domain for the integral of the IP
   * QuadratureFamily quad - quadrature rule to use for the integral of the IP
   * mutable Expr t_ - time at this level is a Sundance::Paramter; mutable allows you to
   *                   affect t_ in functions that are labled as const
   * Expr q_ - Holds the volume force term
   * bool MatrixAndTensorInFile_ - True if A and T are on hand; false if they need to be created
   */
  CellFilter interior_;
  CellFilter boundary_;
  Teuchos::Array<Expr> phi_;
  Expr uB_;
  Expr q_;
  mutable Expr t_;
  double deltat_;
  Mesh mesh_;
  Teuchos::Array<FunctionalEvaluator> forceIP_;
  bool MatrixAndTensorInFile_;
  QuadratureFamily quad_;
  double nu_;
  mutable Expr tNext_;



  
  /********************************************************************************
   * A_IP peforms the IP -(phi_i, uB*(grad*phi_j)) - (phi_i, phi_j*(grad*uB_)) - nu*(grad*f, grad*g)
   ********************************************************************************/
  double A_IP(Expr phi_i, Expr phi_j);

  /********************************************************************************
   * tensorIP peforms the IP -(f, (h*grad)*g)
   ********************************************************************************/
  double tensorIP(Expr f, Expr g, Expr h);

}; // End of MMSQuadODE class

#endif
