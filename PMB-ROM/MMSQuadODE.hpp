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


#include "VientoSnapshotIO.hpp" //For snapshotToMatrix
#include "denseSerialMatrixIO.hpp"
#include "QuadraticODERHSBase.hpp"
#include "MathematicaConverter.hpp"

/** */
class MMSQuadODE : public QuadraticODERHSBase
{
public:
  /**Constructor*/
  MMSQuadODE(Teuchos::Array<Expr> phi, Mesh mesh, bool MatrixAndTensorInFile = false, int verbosity = 1, int quadOrder = 6);

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
  Teuchos::Array<Expr> phi_; 
  Mesh mesh_;
  bool MatrixAndTensorInFile_;
  QuadratureFamily quad_;
  mutable Expr t_;
  Expr q_;

  /********************************************************************************
   * gradIP peforms the IP (grad*f, grad*g)
   ********************************************************************************/
  double gradIP(Expr f, Expr g);

  /********************************************************************************
   * tensorIP peforms the IP (f, (h*grad)*g)
   ********************************************************************************/
  double tensorIP(Expr f, Expr g, Expr h);

}; // End of MMSQuadODE class

#endif
