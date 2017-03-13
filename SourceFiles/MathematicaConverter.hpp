#ifndef MATHEMATICA_CONVERTER_HPP
#define MATHEMATICA_CONVERTER_HPP

#include "Sundance.hpp"

// Things necessary for mathematica to be understood by c++
const double Pi = 4.0*atan(1.0);
Expr Cos(const Expr& x) {return cos(x);}
Expr Sin(const Expr& x) {return sin(x);}
Expr Power(const Expr& x, const double& p) {return pow(x,p);}
Expr Power(const double& x, const Expr& p) {return exp(p*log(x));}
const double E = exp(1.0);
Expr Sqrt(const Expr& x) {return sqrt(x);}
double Sqrt(const double& x) {return sqrt(x);}

#endif
