#ifndef MyNLO_HPP
#define MyNLO_HPP

#include "PlayaSerialVectorType.hpp"
#include "PlayaDenseSerialMatrix.hpp"
#include "PlayaSerialEpetraAdapter.hpp" // Allows me to convert between Epetra and Serial Vectors
#include "QuadraticODERHSBase.hpp" // Need to change if using a different QuadODE
#include "Sundance.hpp"


/*****************************************************************************************
 *
 * The class MyNLO inherits from NonLinearOperatorBase. As such, it must provide
 * implementation for the pure virtual functions:
 *	Vector< double > getInitialGuess() const
 *	LinearOperator< double> computeJacobianAndFunction(Vector<double> &functionValue) const
 *
 * The purpose of MyNLO is work with the function
 *		F(z) = z - uPrev - (h/2.0)*(f(tPrev,uPrev) + f(tNext,z) )
 *
 * Which is the result of the Trapezoid Rule applied to a nonlinear ODE.
 * uNext = z + uPrev
 *
 *****************************************************************************************/
class MyNLO : public NonlinearOperatorBase<double>
{
public:
  //MyNLO(MMSQuadODE f, double h);
  MyNLO(const RCP<QuadraticODERHSBase>& f, double h);
  
  virtual RCP<NonlinearOperatorBase> getRcp() {return rcp(this);}

  Vector<double> getInitialGuess() const {return uPrev_.copy();}

  void set_tPrev(const double t);

  void set_uPrev(const Vector<double> u) {uPrev_ = u.copy();}

protected:
  LinearOperator<double> computeJacobianAndFunction(Vector<double>& functionValue) const;

private:
  //  MMSQuadODE f_;
  RCP<QuadraticODERHSBase> f_;
  Vector<double> uPrev_;
  double h_;
  double tPrev_;
  double tNext_;
};

#endif
