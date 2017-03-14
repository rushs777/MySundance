#ifndef ODERHSBase_HPP
#define ODERHSBase_HPP


#include "PlayaDenseSerialMatrix.hpp"
#include "PlayaSerialVectorType.hpp"
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaVectorSpaceImpl.hpp"
//#include "PlayaSerialEpetraAdapter.hpp" // Allows me to convert between Epetra and Serial Vectors
#include "PlayaOut.hpp"
//#include <vector>


#include "Sundance.hpp"

//using std::vector;
//using std::cout;
using std::endl;
using namespace Teuchos;
using namespace Playa;
using namespace PlayaExprTemplates;

/** Base class for the RHS of a system of ODEs */
class ODERHSBase
{
public:
/**Constructor */
  ODERHSBase(int n);

/**Returns the VectorSpace that the LinearOperator objects are built off*/
  const VectorSpace<double>& space() const;

/**Calculates the RHS of the ODE u_t = f(t,u) */
  virtual Vector<double> eval1(const double& t,
			       const Vector<double>& u,
			       LinearOperator<double>& J) const = 0 ;
private:
/** Holds the VectorSpace that the LinearOperator objects will be built off*/
  VectorSpace<double> vs_;
};


#endif
