#ifndef velocityROM_CLASS_HPP
#define velocityROM_CLASS_HPP

#include "Sundance.hpp"
#include "PlayaSerialEpetraAdapter.hpp"
#include "PlayaSerialVectorType.hpp"
#include "PlayaDenseSerialMatrix.hpp"
#include "denseSerialMatrixIO.hpp"

#include "VientoSnapshotIO.hpp" //For snapshotToMatrix
//#include "PlayaSVD.hpp"
#include "POD_SVD.hpp"
#include "MathematicaConverter.hpp"

//For using newton-armijo
//#include "MMSQuadODE.hpp" //Needs to be changed for further work
#include "NSEProjectedODE.hpp"
#include "MyNLO.hpp"
#include "PlayaDenseLUSolver.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "PlayaNewtonArmijoSolverImpl.hpp"
#include "PlayaNonlinearSolver.hpp"
#include "PlayaNonlinearSolverBuilder.hpp"

/**The class velocityROM generates the reduced order velocity using velocity POD modes*/
class velocityROM
{
public:
  /**Constructor*/
  velocityROM(const string POD_base_dir, string nlParamFile, const DiscreteSpace& ds,
	      Expr u0, Expr forceTerm, Expr t, int nSteps, double deltat, double nu,
	      double tolerance = .999, int verbosity = 1);

  /**initialize will set up the fluctuation matrix Wprime and the fluctuation velocity POD basis functions*/
  void initialize(const string snapshotFilePrefix);

  /**get_phi will return the fluctuation velocity POD basis functions*/
  Array<Expr> get_phi() {return phi_;}

  /**get_uB will return the time mean for the velocity*/
  Expr get_uB() {return uB_;}

  /**generate_alpha will calculate alpha(t_m) for m=0:nSteps*/
  void generate_alpha();

  /**Return the Array containing the alpha(t_m) values */
  Array<Vector<double> > get_alpha() {return alpha_;}

  /**Make the reduced-order pressure term from ubar_, alpha_, and phi_*/
  Array<Expr> get_uRO();

  /** Writes the time-dependent vector from the ODE constraint to filename */
  void write_b(const string filename);

  /** Writes the time-dependent coefficient function alpha to filename */
  void write_alpha(const string filename);

  /** alphaErrorCheck calculates the exact projected form of the time-dependent functions
   * and returns a Vector<double> object which contains the 2norm of the error 
   * between alphaEx and alpha_ at each time step;
   * i.e. alphaError[i] = || alphaEx(t_i) - alpha_(t_i)|| 
   * IMPORTANT: The Expr uExact has to have the same Sundance::Parameter for time
   * as the one that was passed into the constructor.
   */
  Vector<double> alphaErrorCheck(Expr uExact);

  /** */
  Vector<double> velocityErrorCheck(Expr uExact);


private:
  string POD_base_dir_;
  string nlParamFile_;
  DiscreteSpace ds_;
  Expr u0_;
  Expr forceTerm_;
  Expr t_;
  int nSteps_;
  double deltat_;
  double nu_;
  double tol_;
  int verbosity_;
  int R_;
  Expr uB_;
  //LinearOperator<double> Wprime_;
  Array<Expr> phi_;
  Array<Vector<double> > alpha_;


};

#endif
