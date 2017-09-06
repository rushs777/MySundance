#ifndef ODECO_CLASS_HPP
#define ODECO_CLASS_HPP

// Sundance Includes
#include "Sundance.hpp"
#include "PlayaDenseSerialMatrix.hpp"
#include "VientoSnapshotIO.hpp" // for readSnap
#include "denseSerialMatrixIO.hpp"

// Local Includes
#include "MathematicaConverter.hpp"


/**
 * ODECO is a class designed to perform ODE-Constrained Optimization
 */
class ODECO
{
public:
  /** Constructor */
  ODECO(string ROM_base_dir, string matrixAndTensorDir, Mesh spatialMesh, int verbosity=0);

  void initialize();

  Expr get_b();

  Expr get_phi();

  Expr get_uB();

  Expr get_A();

  Expr get_At();

  Expr get_T();


  double errorCheck(Expr alphaOPT);
  
private:
  string ROM_base_dir_;
  string matrixAndTensorDir_;
  Mesh spatialMesh_;
  Mesh timeMesh_;
  VectorType<double> epetraVecType_;
  int Ru_;
  BasisFamily time_basis_;
  Expr b_;
  Array<Expr> phi_;
  Expr uB_;
  Expr A_;
  Expr At_;
  Expr T_;
  int verbosity_;

  void initialize_b();

  void initialize_phi();

  void initialize_uB();

  void initialize_A_and_T();
};



#endif
