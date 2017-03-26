#ifndef velocityROM_CLASS_HPP
#define velocityROM_CLASS_HPP

#include "Sundance.hpp"
#include "PlayaSerialEpetraAdapter.hpp"
#include "PlayaSerialVectorType.hpp"
#include "PlayaDenseSerialMatrix.hpp"
#include "PlayaDenseLUSolver.hpp"

#include "VientoSnapshotIO.hpp" //For snapshotToMatrix
#include "PlayaSVD.hpp"
#include "MathematicaConverter.hpp"

//For using newton-armijo
#include "MMSQuadODE.hpp" //Needs to be changed for further work
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
  velocityROM(const string snapshotFilename, string nlParamFile, const DiscreteSpace& ds, Expr u0, Expr forceTerm, Expr t, int nSteps, double deltat, double tolerance = .999, int verbosity = 1);

  /**initialize will set up the fluctuation matrix Wprime and the fluctuation velocity POD basis functions*/
  void initialize();

  /**get_phi will return the fluctuation velocity POD basis functions*/
  Array<Expr> get_phi() {return phi_;}

  /**generate_alpha will calculate alpha(t_m) for m=0:nSteps*/
  void generate_alpha();

  /**Return the Array containing the alpha(t_m) values */
  Array<Vector<double> > get_alpha() {return alpha_;}

  /**Make the reduced-order pressure term from ubar_, alpha_, and phi_*/
  Array<Expr> get_uRO();

private:
  string snapshotFilename_;
  string nlParamFile_;
  DiscreteSpace ds_;
  Expr u0_;
  Expr forceTerm_;
  Expr t_;
  int nSteps_;
  double deltat_;
  double tol_;
  int verbosity_;
  int R_;
  Expr ubar_;
  LinearOperator<double> Wprime_;
  Array<Expr> phi_;
  Array<Vector<double> > alpha_;


};

#endif
