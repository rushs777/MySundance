#ifndef QuadraticODERHSBase_HPP
#define QuadraticODERHSBase_HPP

#include "PlayaSerialVectorType.hpp"
#include "PlayaDenseSerialMatrix.hpp"
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaVectorSpaceImpl.hpp"
//#include "PlayaSerialEpetraAdapter.hpp" // Allows me to convert between Epetra and Serial Vectors
#include "PlayaOut.hpp"
//#include <vector>

#include "ODERHSBase.hpp"


#include "Sundance.hpp"

//using std::vector;
//using std::cout;
using std::endl;
using namespace Teuchos;
using namespace Playa;
using namespace PlayaExprTemplates;


/** */
class QuadraticODERHSBase : public ODERHSBase
{
public:
/**Constructor*/
  QuadraticODERHSBase(int n, int verbosity = 1);

/**Initializes the class with it's required matrix and tensor*/
  void initialize();

/**Evaluates the force term in the ODE*/
  virtual Vector<double> evalForceTerm(const double& t) const = 0 ;

/**Returns the RHS of the ODE u_t = f(t,u); calculates the Jacobian of f*/
  Vector<double> eval1(const double& t,
		       const Vector<double>& u,
		       LinearOperator<double>& J) const;

/**Pure virtual function for the problem specific task of calculating the values 
*  for the matrix A and the tensor T
*/
  virtual void fillMatrixAndTensor(RCP<DenseSerialMatrix>& A,
				   Array<RCP<DenseSerialMatrix> >& T) = 0 ;

/** Returns the verbosity level for this RHS */
  int getVerbosity() const ; // KRL: this should be const



private:
/**Jacobain of f*/
  mutable LinearOperator<double> J_;
/** */
  LinearOperator<double> A_;
/** */
  Array<LinearOperator<double> > T_;
/** */
  bool initialized_;
/** */
  int verbosity_;
};



#endif
