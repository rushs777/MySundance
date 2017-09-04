#ifndef KKT_BASE_HPP
#define KKT_BASE_HPP

// Sundance Includes
#include "Sundance.hpp"
#include "PlayaDenseSerialMatrix.hpp"
#include "VientoSnapshotIO.hpp" // for readSnap
#include "denseSerialMatrixIO.hpp"

#include "ODECO.hpp"

// Local Includes
#include "MathematicaConverter.hpp"


/**
 * KKTBase assumes the time_basis is Lagrange(1)
 */
class KKTBase
{
public:
  /** Constructor */
  KKTBase(Mesh timeMesh, int Ru);

  /** Establish the state equation and it's boundary condition. 
      Must be done by hand for each problem type */
  virtual void stateEqn() = 0;

  /** Establish the adjoint equation and it's boundary condition 
      Must be done by hand for each problem type */
  virtual void adjointEqn() = 0;

  /** Establish the design equation and it's boundary condition 
      Must be done by hand for each problem type */
  virtual void designEqn() = 0;

  /** 
   * Solve the KKT system .
   * Returns the array of solved values as an Expr
   */
  virtual Expr solve() = 0;


protected:
  /** 
   * State Variable: alpha_
   * Adjoint Variable: lambda_
   * Design Variable: p_
   * Dimension of each variable: Ru_
   */
  Expr alpha_;
  Expr lambda_;
  Expr p_;
  Expr alphaHat_;
  Expr lambdaHat_;
  Expr pHat_;
  int Ru_;
  DiscreteSpace ODECO_DS_;
  Expr dt_;
  Expr stateEqn_;
  Expr stateBC_;
  Expr adjointEqn_;
  Expr adjointBC_;
  Expr designEqn_;
  Expr designBC_;
  CellFilter interior_;
};



#endif
