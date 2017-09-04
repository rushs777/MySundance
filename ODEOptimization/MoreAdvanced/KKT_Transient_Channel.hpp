#ifndef KKT_TRANSIENT_CHANNEL_HPP
#define KKT_TRANSIENT_CHANNEL_HPP

#include "KKTBase.hpp"

/**
 * This derived class models the KKT system arising from a Transient Channel Problem
 */
class KKT_Transient_Channel : public KKTBase
{
public:
  /** Constructor */
  KKT_Transient_Channel(Mesh timeMesh, int Ru);

  /** Establish the state equation and it's boundary condition */
  virtual void stateEqn();

  /** Establish the adjoint equation and it's boundary condition */
  virtual void adjointEqn();

  /** Establish the design equation and it's boundary condition */
  virtual void designEqn();

  /** 
   * Solve the KKT system .
   * Returns the array of solved values as an Expr
   */
  virtual Expr solve();
  


private:

  
};

#endif
