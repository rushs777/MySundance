#ifndef Playa_SVD_HPP
#define Playa_SVD_HPP

// Code
#include "PlayaLinearOperatorDecl.hpp"// For LinearOperator
#include "Sundance.hpp"

namespace Playa
{

  void POD(const LinearOperator<double> &W, Vector<double> &lambda, LinearOperator<double> &Alpha, LinearOperator<double> &Phi, Sundance::DiscreteSpace &ds, int debug = 1);
}
#endif
