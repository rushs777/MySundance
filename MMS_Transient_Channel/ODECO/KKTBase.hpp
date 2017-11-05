#ifndef KKT_BASE_HPP
#define KKT_BASE_HPP

// Sundance Includes
#include "Sundance.hpp"
#include "PlayaDenseSerialMatrix.hpp"
#include "VientoSnapshotIO.hpp" // for readSnap
#include "denseSerialMatrixIO.hpp"
#include "POD_SVD.hpp"


// Local Includes
//#include "MathematicaConverter.hpp"


/**
 * KKTBase is a base class for solving the KKT system of equations.
 * The state, adjoint, and design equations are specified in derived classes.
 * It assumes the time basis is Lagrange(1).
 */
class KKTBase
{
public:
  /** Constructor */
  KKTBase(string POD_DataDir, Mesh spatialMesh, Mesh timeMesh, double tFinal, int nSteps, int verbosity=0);

  /** Initialize will read in from file the values for necessary matrices, tensor, and vector-valued functions */
  void initialize();

  /** 
   * Solve the KKT system .
   * Returns the array of solved values as an Expr
   */
  virtual Expr solve(string solverXML, double eta_design, double eta_reg=0) = 0;

  /**
   * Return the value of alphaOPT_
   */
  Expr get_alpha() {return alphaOPT_;}

  /**
   * Return the value of uOPT_
   */
  Array<Expr> get_uOPT() {return uOPT_;}

protected:
  /** 
   * State Variable: alpha_
   * Adjoint Variable: lambda_
   * Design Variable: p_
   * Dimension of each variable: Ru_
   * The test function form of the variables are fooHat
   */
  // Initialized in the constructor
  string POD_DataDir_;
  Mesh spatialMesh_;
  Mesh timeMesh_;
  double tFinal_;
  int nSteps_;
  int verbosity_;
  Expr dt_;
  CellFilter interior_;

  // Initialized in initialize()
  VectorType<double> epetraVecType_;

  // Initialized in initialize_b()
  BasisFamily time_basis_;
  Expr b_;
  //  DiscreteSpace time_DS_; Currently not used anywhere else, and so not a class variable

  // Initialized in initialize_phi()
  int Ru_;
  Array<Expr> phi_;

  // Initialized in initialize_uB()
  Expr uB_;

  // Initialized in initialize_A_and_T()
  Expr A_;
  Expr At_;
  Expr T_;
  
  // Initialized in initialize_vars()
  Expr alpha_;
  Expr lambda_;
  Expr p_;
  Expr alphaHat_;
  Expr lambdaHat_;
  Expr pHat_;

  // Initialized in derived classes
  Expr stateEqn_;
  Expr stateBC_;
  Expr adjointEqn_;
  Expr adjointBC_;
  Expr designEqn_;
  Expr designBC_;
  Expr alphaOPT_;
  Array<Expr> uOPT_;



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
