#ifndef KKT_TRANSIENT_CHANNEL_HPP
#define KKT_TRANSIENT_CHANNEL_HPP

#include "KKTBase.hpp"
#include "sensorData.hpp"

/**
 * This derived class models the KKT system arising from a Transient Channel Problem
 */
class KKT_Transient_Channel : public KKTBase
{
public:
  /** Constructor */
  KKT_Transient_Channel(string ROM_base_dir, Mesh spatialMesh, Mesh timeMesh, sensorData dataClass, double tFinal, int quadOrder, int verbosity=0);



  /** 
   * Solve the KKT system .
   * Returns the array of solved values as an Expr
   */
  Expr solve(string solverXML, double eta_design, double eta_reg);

  /** 
   * In the instance where the parameter space only has one seleciton of values,
   * compare alphaOPT to alphaROM
   */
  double errorCheck(Expr alphaOPT);
  


private:
  sensorData dataClass_;
  QuadratureFamily quad_;


  /** Set up the regularization term */
  virtual void regularizationTerm(double eta);

  /** Establish the state equation and it's boundary condition */
  virtual void stateEqn();

  /** Establish the adjoint equation and it's boundary condition */
  virtual void adjointEqn();

  /** Establish the design equation and it's boundary condition */
  virtual void designEqn(double eta);

  
};

#endif
