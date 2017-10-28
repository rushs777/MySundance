#ifndef NSEProjectedODE_2D_HPP
#define NSEProjectedODE_2D_HPP

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
#include "meshAddOns.hpp"
#include "MathematicaConverter.hpp"

/** NSEProjectedODE creates the necessary components from projecting the NSE onto the space spanned by a set of basis functions */
class NSEProjectedODE : public QuadraticODERHSBase
{
public:
  /**Constructor*/
  NSEProjectedODE(Teuchos::Array<Expr> phi, Expr uB, Expr q, Expr t, double deltat, Mesh mesh, bool MatrixAndTensorInFile = false, int verbosity = 1, int quadOrder = 6);

  /**Project the force term onto the basis*/
  Vector<double> evalForceTerm(const double& t) const;

  /**Create A and T according to this problem*/
  void fillMatrixAndTensor(RCP<DenseSerialMatrix>& A, Array<RCP<DenseSerialMatrix> >& T);

  /**Write the DiscreteFunctionVectors for the time-dependent function b to file*/
  void write_b(int nSteps);

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
  string fileDir_;

  /** Avoid expensive calculations of q by discretizing q onto a discrete space, once
   * per timestep. */
  mutable map<double, Expr> qCache_;
  mutable Expr qDiscrete_;
  L2Projector qProj_;
  void updateQ(const double& t) const ;



  
  /********************************************************************************
   * A_IP peforms the IP -(phi_i, uB*(grad*phi_j)) - (phi_i, phi_j*(grad*uB_)) - nu*(grad*f, grad*g)
   ********************************************************************************/
  double A_IP(Expr phi_i, Expr phi_j);

  /********************************************************************************
   * tensorIP peforms the IP -(f, (h*grad)*g)
   ********************************************************************************/
  double tensorIP(Expr f, Expr g, Expr h);

}; // End of MMSQuadODE class


/**
 * This class is designed to give a CELL_PREDICATE test that
 * returns true for all Points that do not have an x-coordinate 
 * equal to dim0_max_
 */
class outflowEdgeTest: public CellPredicateFunctorBase,
		       public Playa::Handleable<CellPredicateFunctorBase>
{
public:
  outflowEdgeTest(const double& dim0_max) : CellPredicateFunctorBase("outflowEdgeTest"), dim0_max_(dim0_max) {}
  virtual ~outflowEdgeTest() {}
  virtual bool operator()(const Point& x) const {return fabs(x[0] - dim0_max_) > 1.0e-10;}
  GET_RCP(CellPredicateFunctorBase);

private:
  double dim0_max_;
};




#endif
