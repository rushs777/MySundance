#ifndef MATHEMATICA_CONVERTER_HPP
#define MATHEMATICA_CONVERTER_HPP

#include "Sundance.hpp"

// Things necessary for mathematica to be understood by c++
const double Pi = 4.0*atan(1.0);
inline Expr Cos(const Expr& x) {return cos(x);} 
inline Expr Sin(const Expr& x) {return sin(x);}
inline Expr Csc(const Expr& x) {return (1.0/sin(x));}
inline Expr Power(const Expr& x, const double& p) {return pow(x,p);}
inline Expr Power(const double& x, const Expr& p) {return exp(p*log(x));}
const double E = exp(1.0);
inline Expr Sqrt(const Expr& x) {return sqrt(x);}
inline double Sqrt(const double& x) {return sqrt(x);}

// KRL: declare these functions inline to avoid duplicate symbol errors at link time

#endif
