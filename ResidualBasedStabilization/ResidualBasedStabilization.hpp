#ifndef ResidualBasedStabilization_CLASS_HPP
#define ResidualBasedStabilization_CLASS_HPP

#include "Sundance.hpp"
#include "PlayaSerialEpetraAdapter.hpp"
#include "PlayaSerialVectorType.hpp"
#include "PlayaDenseSerialMatrix.hpp"

// Local Files
#include "PlayaSVD.hpp"
#include "MathematicaConverter.hpp"
#include "velocityROM.hpp"
#include "PMB.hpp"

/**The class ResidualBasedStabilization attempts to improve the ROM by adding 
 * in POD modes from the residual of the NSE
 */
class ResidualBasedStabilization
{
public:
  /**Constructor*/
  ResidualBasedStabilization(velocityROM velROM, PMB pressROM, Expr forceTerm, Expr t, double deltat, int K, int verbosity=1);

  /**initialize will calcluate the residuals for the momentum and continuity equations (Rm, Rc)
   * and organize their results into matrices for which we perform a POD
   */
  void initialize();

  /** Generate the orthonormal basis functions for both stabilized velocity and pressure*/
  void generateBasisFunctions();

  /** Get the stabilized reduced-order velocity*/
  Array<Expr> get_uRO();

  /** Get the stabilized reduced-order pressure*/
  Array<Expr> get_pRO();

private:
  /** Calculate the residual for the momentum equation */
  Vector<double> calculateRm(int timeIndex);

  /** Calculate the residual for the continuity equation */
  Vector<double> calculateRc(int timeIndex);

  /** Perform Gram-Schmidt on the basis functions to get an orthonormal set */
  Array<Expr> GramSchmidt(Array<Expr> v, Mesh mesh, int quadOrder);

  velocityROM velROM_;
  PMB pressROM_;
  Expr forceTerm_;
  Expr t_;
  double deltat_;
  int K_;
  int verbosity_;
  Array<Expr> uRO_;
  Array<Expr> pRO_;
  DiscreteSpace Rmds_;
  DiscreteSpace Rcds_;
  double nu_;
  Array<Expr> zeta_;
  Array<Expr> eta;
  Array<Expr> velBasisFuncs_;
  Array<Expr> pressBasisFuncs_;



  
}






#endif
