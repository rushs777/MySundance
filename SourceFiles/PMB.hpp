#ifndef PMB_CLASS_HPP
#define PMB_CLASS_HPP

#include "Sundance.hpp"
#include "PlayaSerialEpetraAdapter.hpp"
#include "PlayaSerialVectorType.hpp"
#include "PlayaDenseSerialMatrix.hpp"
#include "PlayaDenseLUSolver.hpp"

#include "VientoSnapshotIO.hpp" //For snapshotToMatrix
#include "PlayaSVD.hpp"

/**The class PMB generates the pressure based off of pressure POD modes*/
class PMB
{
public:
  /**Constructor*/
  PMB(const string filename, const Array<Expr>& uRO, Expr forceTerm, Expr t, const DiscreteSpace& ds, double deltat, double tolerance = .999, int verbosity = 1);

  /**initialize will set up the fluctuation matrix Qprime and the fluctuation pressure POD basis functions*/
  void initialize();

  /**get_psi will return the pressure POD basis functions*/
  Array<Expr> get_psi() {return psi_;}

  /**generate_beta will calculate beta(t_m) for m=0:nSteps*/
  void generate_beta();

  /**Return the Array containing the beta(t_m) values */
  Array<Vector<double> > get_beta() {return beta_;}

  /**Make the reduced-order pressure term from pbar_ and beta*/
  Array<Expr> get_pRO();

private:
  string filename_;
  Array<Expr> uRO_;
  Expr forceTerm_;
  Expr t_;
  DiscreteSpace ds_;
  double deltat_;
  int M_;
  double tol_;
  int verbosity_;
  int R_;
  Expr pbar_;
  LinearOperator<double> Qprime_;
  Array<Expr> psi_;
  Array<Vector<double> > beta_;


  /********************************************************************************
   * gradIP peforms the IP (grad*f, grad*g)
   ********************************************************************************/
  double gradIP(Expr f, Expr g);

  double L2IP(Expr f, Expr g);

};

#endif
