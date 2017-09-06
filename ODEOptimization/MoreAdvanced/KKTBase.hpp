#ifndef KKT_BASE_HPP
#define KKT_BASE_HPP

// Sundance Includes
#include "Sundance.hpp"
#include "PlayaDenseSerialMatrix.hpp"
#include "VientoSnapshotIO.hpp" // for readSnap
#include "denseSerialMatrixIO.hpp"


// Local Includes
#include "MathematicaConverter.hpp"


/**
 * KKTBase is a base class for solving the KKT system of equations.
 * It assumes the time basis is Lagrange(1).
 */
class KKTBase
{
public:
  /** Constructor */
  KKTBase(string ROM_base_dir, Mesh spatialMesh, Mesh timeMesh, double tFinal, int verbosity=0);

  /** Initialize will read in from file the values for necessary matrices, tensor, and vector-valued functions */
  void initialize();

  /** 
   * Solve the KKT system .
   * Returns the array of solved values as an Expr
   */
  virtual Expr solve(string solverXML, double eta_design, double eta_reg=0) = 0;


protected:
  /** 
   * State Variable: alpha_
   * Adjoint Variable: lambda_
   * Design Variable: p_
   * Dimension of each variable: Ru_
   * The test function form of the variables are fooHat
   */
  Expr alpha_;
  Expr lambda_;
  Expr p_;
  Expr alphaHat_;
  Expr lambdaHat_;
  Expr pHat_;
  int Ru_;
  Expr dt_;
  Expr stateEqn_;
  Expr stateBC_;
  Expr adjointEqn_;
  Expr adjointBC_;
  Expr designEqn_;
  Expr designBC_;
  CellFilter interior_;


  string ROM_base_dir_;
  Mesh spatialMesh_;
  Mesh timeMesh_;
  double tFinal_;
  VectorType<double> epetraVecType_;
  BasisFamily time_basis_;
  //  DiscreteSpace time_DS_;
  Expr b_;
  Array<Expr> phi_;
  Expr uB_;
  Expr A_;
  Expr At_;
  Expr T_;
  int verbosity_;

  /** Establish the state equation and it's boundary condition. 
      Must be done by hand for each problem type */
  virtual void stateEqn() = 0;

  /** Establish the adjoint equation and it's boundary condition 
      Must be done by hand for each problem type */
  virtual void adjointEqn() = 0;

  /** Establish the design equation and it's boundary condition 
      Must be done by hand for each problem type */
  virtual void designEqn(double eta) = 0;

  void initialize_b();

  void initialize_phi();

  void initialize_uB();

  void initialize_A_and_T();

  void initialize_vars();

};



#endif
