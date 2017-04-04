#include "VientoPCDStepper.hpp"
#include "VientoErrorChecker.hpp"
#include "Sundance.hpp"
#include "VientoDefaultOutputManager.hpp"



using namespace Sundance;
using namespace Viento;
using Sundance::List;

// Things necessary for mathematica generated fs to be understood
#include "MathematicaConverter.hpp"

/** 
 * MMS Inflow with a projection method
 */

CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-2.0*M_PI) < 1.0e-10;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-M_PI) < 1.0e-10;})

CELL_PREDICATE(CornerPointTest, {return fabs(x[1]) < 1.0e-10 && fabs(x[0])<1.0e-10;})

Expr uExFunc(const Expr& t)
{
  Expr x = new CoordExpr(0);
  Expr y = new CoordExpr(1);

  return List(1 + (x*(-(Cos(0.5 - t)*Cos(y)*
             Sin(0.5 - t)*Sin(x))/5.\
          - (2*Cos(2*(-0.5 + t))*
            Cos(2*y)*Sin(0.5 - t)*
            Sin(x))/5. + 
         (Cos(0.5 - t)*Cos(y)*
            Sin(2*(-0.5 + t))*
            Sin(2*x))/3. + 
         (2*Cos(2*(-0.5 + t))*
            Cos(2*y)*
            Sin(2*(-0.5 + t))*
            Sin(2*x))/5. + 
         (Cos(0.5 - t)*Cos(y)*
            Sin(3*(-0.5 + t))*
            Sin(3*x))/6. + 
         (2*Cos(2*(-0.5 + t))*
            Cos(2*y)*
            Sin(3*(-0.5 + t))*
            Sin(3*x))/7.))/(1 + x),
   -((x*(-(Cos(0.5 - t)*Cos(x)*
              Sin(0.5 - t)*Sin(y))/5.
             + (2*Cos(0.5 - t)*
              Cos(2*x)*
              Sin(2*(-0.5 + t))*
              Sin(y))/3. + 
           (Cos(0.5 - t)*Cos(3*x)*
              Sin(3*(-0.5 + t))*
              Sin(y))/2. - 
           (Cos(2*(-0.5 + t))*Cos(x)*
              Sin(0.5 - t)*Sin(2*y))/
            5. + 
           (2*Cos(2*(-0.5 + t))*
              Cos(2*x)*
              Sin(2*(-0.5 + t))*
              Sin(2*y))/5. + 
           (3*Cos(2*(-0.5 + t))*
              Cos(3*x)*
              Sin(3*(-0.5 + t))*
              Sin(2*y))/7.))/(1 + x))
      + (x*(-(Cos(0.5 - t)*
             Sin(0.5 - t)*Sin(x)*
             Sin(y))/5. + 
         (Cos(0.5 - t)*
            Sin(2*(-0.5 + t))*
            Sin(2*x)*Sin(y))/3. + 
         (Cos(0.5 - t)*
            Sin(3*(-0.5 + t))*
            Sin(3*x)*Sin(y))/6. - 
         (Cos(2*(-0.5 + t))*
            Sin(0.5 - t)*Sin(x)*
            Sin(2*y))/5. + 
         (Cos(2*(-0.5 + t))*
            Sin(2*(-0.5 + t))*
            Sin(2*x)*Sin(2*y))/5. + 
         (Cos(2*(-0.5 + t))*
            Sin(3*(-0.5 + t))*
            Sin(3*x)*Sin(2*y))/7.))/
     Power(1 + x,2) - 
    (-(Cos(0.5 - t)*Sin(0.5 - t)*
           Sin(x)*Sin(y))/5. + 
       (Cos(0.5 - t)*
          Sin(2*(-0.5 + t))*Sin(2*x)*
          Sin(y))/3. + 
       (Cos(0.5 - t)*
          Sin(3*(-0.5 + t))*Sin(3*x)*
          Sin(y))/6. - 
       (Cos(2*(-0.5 + t))*
          Sin(0.5 - t)*Sin(x)*
          Sin(2*y))/5. + 
       (Cos(2*(-0.5 + t))*
          Sin(2*(-0.5 + t))*Sin(2*x)*
          Sin(2*y))/5. + 
       (Cos(2*(-0.5 + t))*
          Sin(3*(-0.5 + t))*Sin(3*x)*
	Sin(2*y))/7.)/(1 + x));
}

Expr pExFunc(const Expr& t)
{
  Expr x = new CoordExpr(0);
  Expr y = new CoordExpr(1);

  return 2*Pi - x + Cos(y)*Sin(x);
}

Expr qExFunc(const Expr& t)
{
  Expr x = new CoordExpr(0);
  Expr y = new CoordExpr(1);
  double nu = 1.0;

  return List((-Power(1 + x,3) + 
      Power(1 + x,3)*Cos(x)*Cos(y) + 
      ((7*x*Cos(0.5 - t)*Cos(y)*
            (-6*Sin(0.5 - t)*
             Sin(x) + 
              10*Sin(2*(-0.5 + t))*
              Sin(2*x) + 
              5*Sin(3*(-0.5 + t))*
              Sin(3*x)) + 
           6*(7*
              (5 + 5*x + 
              x*Cos(2*y)*
              Sin(4*(-0.5 + t))*
              Sin(2*x)) - 
              2*x*Cos(2*(-0.5 + t))*
              Cos(2*y)*
              (7*Sin(0.5 - t)*
             Sin(x) - 
              5*Sin(3*(-0.5 + t))*
              Sin(3*x))))*
         (7*Cos(0.5 - t)*Cos(y)*
            (-6*x*(1 + x)*Cos(x)*
              Sin(0.5 - t) + 
              20*x*(1 + x)*Cos(2*x)*
              Sin(2*(-0.5 + t)) + 
              15*x*Cos(3*x)*
              Sin(3*(-0.5 + t)) + 
              15*Power(x,2)*Cos(3*x)*
              Sin(3*(-0.5 + t)) - 
              6*Sin(0.5 - t)*
             Sin(x) + 
              10*Sin(2*(-0.5 + t))*
              Sin(2*x) + 
              5*Sin(3*(-0.5 + t))*
              Sin(3*x)) + 
           6*Cos(2*y)*
            (7*Sin(4*(-0.5 + t))*
              (2*x*(1 + x)*
             Cos(2*x) + Sin(2*x)) + 
              2*Cos(2*(-0.5 + t))*
              (-7*x*(1 + x)*Cos(x)*
              Sin(0.5 - t) + 
              15*x*(1 + x)*Cos(3*x)*
              Sin(3*(-0.5 + t)) - 
              7*Sin(0.5 - t)*
             Sin(x) + 
              5*Sin(3*(-0.5 + t))*
              Sin(3*x)))))/44100. + 
      (x*Power(1 + x,2)*
         (42*Power(Cos(0.5 - t),2)*
            Cos(y)*Sin(x) + 
           7*Cos(y)*Sin(0.5 - t)*
            (-6*Sin(0.5 - t) + 
              5*
              (4*Cos(x)*
              Sin(2*(-0.5 + t)) + 
              (1 + 2*Cos(2*x))*
              Sin(3*(-0.5 + t))))*
            Sin(x) + 
           7*Cos(0.5 - t)*
            (4*Cos(2*(-0.5 + t))*
              (3*Cos(2*y)*Sin(x) + 
              5*Cos(y)*Sin(2*x)) + 
              15*Cos(3*(-0.5 + t))*
              Cos(y)*Sin(3*x)) + 
           12*Cos(2*y)*
            (14*Sin(0.5 - t)*
              Sin(2*(-0.5 + t))*
              Sin(x) + 
              14*
              Power(Cos(2*
             (-0.5 + t)),2)*Sin(2*x)\
              + 15*Cos(2*(-0.5 + t))*
              Cos(3*(-0.5 + t))*
              Sin(3*x) - 
              2*Sin(2*(-0.5 + t))*
              (7*Sin(2*(-0.5 + t))*
              Sin(2*x) + 
              5*Sin(3*(-0.5 + t))*
              Sin(3*x)))))/210. + 
      (nu*(7*Cos(0.5 - t)*Cos(y)*
            (6*(1 + x)*Cos(x)*
              Sin(0.5 - t) - 
              20*(1 + x)*Cos(2*x)*
              Sin(2*(-0.5 + t)) - 
              15*Cos(3*x)*
              Sin(3*(-0.5 + t)) - 
              15*x*Cos(3*x)*
              Sin(3*(-0.5 + t)) - 
              6*Sin(0.5 - t)*
             Sin(x) - 
              6*x*Sin(0.5 - t)*
              Sin(x) - 
              6*Power(x,3)*
              Sin(0.5 - t)*Sin(x) + 
              10*Sin(2*(-0.5 + t))*
              Sin(2*x) + 
              25*x*Sin(2*(-0.5 + t))*
              Sin(2*x) + 
              50*Power(x,2)*
              Sin(2*(-0.5 + t))*
              Sin(2*x) + 
              25*Power(x,3)*
              Sin(2*(-0.5 + t))*
              Sin(2*x) + 
              5*Sin(3*(-0.5 + t))*
              Sin(3*x) + 
              25*x*Sin(3*(-0.5 + t))*
              Sin(3*x) + 
              50*Power(x,2)*
              Sin(3*(-0.5 + t))*
              Sin(3*x) + 
              25*Power(x,3)*
              Sin(3*(-0.5 + t))*
              Sin(3*x)) + 
           6*(7*
              (-(Power(x,2)*Cos(y)*
              Sin(1. - 2.*t)*Sin(x))\
              + (1 + 4*x + 
              8*Power(x,2) + 
              4*Power(x,3))*Cos(2*y)*
              Sin(4*(-0.5 + t))*
              Sin(2*x)) + 
              Cos(2*(-0.5 + t))*
              Cos(2*y)*
              (14*(1 + x)*Cos(x)*
              Sin(0.5 - t) - 
              28*(1 + x)*Cos(2*x)*
              Sin(2*(-0.5 + t)) - 
              30*Cos(3*x)*
              Sin(3*(-0.5 + t)) - 
              30*x*Cos(3*x)*
              Sin(3*(-0.5 + t)) - 
              14*Sin(0.5 - t)*
             Sin(x) - 
              35*x*Sin(0.5 - t)*
              Sin(x) - 
              70*Power(x,2)*
              Sin(0.5 - t)*Sin(x) - 
              35*Power(x,3)*
              Sin(0.5 - t)*Sin(x) + 
              10*Sin(3*(-0.5 + t))*
              Sin(3*x) + 
              65*x*Sin(3*(-0.5 + t))*
              Sin(3*x) + 
              130*Power(x,2)*
              Sin(3*(-0.5 + t))*
              Sin(3*x) + 
              65*Power(x,3)*
              Sin(3*(-0.5 + t))*
              Sin(3*x)))))/105. - 
      (x*(21*Sin(1. - 2.*t)*Sin(x)*
            Sin(y) - 
           35*Cos(0.5 - t)*
            (2*Sin(2*(-0.5 + t))*
              Sin(2*x) + 
              Sin(3*(-0.5 + t))*
              Sin(3*x))*Sin(y) + 
           12*(-7*Sin(4*(-0.5 + t))*
              Sin(2*x) + 
              2*Cos(2*(-0.5 + t))*
              (7*Sin(0.5 - t)*
             Sin(x) - 
              5*Sin(3*(-0.5 + t))*
              Sin(3*x)))*Sin(2*y))*
         (-21*Sin(1. - 2.*t)*Sin(x)*
            Sin(y) + 
           7*Cos(0.5 - t)*
            (-6*x*(1 + x)*Cos(x)*
              Sin(0.5 - t) + 
              5*
              (4*x*(1 + x)*Cos(2*x)*
              Sin(2*(-0.5 + t)) + 
              3*x*(1 + x)*Cos(3*x)*
              Sin(3*(-0.5 + t)) + 
              2*Sin(2*(-0.5 + t))*
              Sin(2*x) + 
              Sin(3*(-0.5 + t))*
              Sin(3*x)))*Sin(y) + 
           3*(7*Sin(4*(-0.5 + t))*
              (2*x*(1 + x)*
             Cos(2*x) + Sin(2*x)) + 
              2*Cos(2*(-0.5 + t))*
              (-7*x*(1 + x)*Cos(x)*
              Sin(0.5 - t) + 
              15*x*(1 + x)*Cos(3*x)*
              Sin(3*(-0.5 + t)) - 
              7*Sin(0.5 - t)*
             Sin(x) + 
              5*Sin(3*(-0.5 + t))*
              Sin(3*x)))*Sin(2*y)))/
       44100.)/Power(1 + x,3),
   (-44100*Power(1 + x,4)*Sin(x)*
       Sin(y) + 
      (-7*Cos(0.5 - t)*Cos(y)*
          (-6*x*(1 + x)*Cos(x)*
             Sin(0.5 - t) + 
            20*x*(1 + x)*Cos(2*x)*
             Sin(2*(-0.5 + t)) + 
            15*x*Cos(3*x)*
             Sin(3*(-0.5 + t)) + 
            15*Power(x,2)*Cos(3*x)*
             Sin(3*(-0.5 + t)) - 
            6*Sin(0.5 - t)*Sin(x) + 
            10*Sin(2*(-0.5 + t))*
             Sin(2*x) + 
            5*Sin(3*(-0.5 + t))*
             Sin(3*x)) - 
         6*Cos(2*y)*
          (7*Sin(4*(-0.5 + t))*
             (2*x*(1 + x)*Cos(2*x) + 
              Sin(2*x)) + 
            2*Cos(2*(-0.5 + t))*
             (-7*x*(1 + x)*Cos(x)*
              Sin(0.5 - t) + 
              15*x*(1 + x)*Cos(3*x)*
              Sin(3*(-0.5 + t)) - 
              7*Sin(0.5 - t)*
             Sin(x) + 
              5*Sin(3*(-0.5 + t))*
              Sin(3*x))))*
       (-7*Cos(0.5 - t)*
          (-6*x*(1 + x)*Cos(x)*
             Sin(0.5 - t) + 
            5*(4*x*(1 + x)*Cos(2*x)*
              Sin(2*(-0.5 + t)) + 
              3*x*(1 + x)*Cos(3*x)*
              Sin(3*(-0.5 + t)) + 
              2*Sin(2*(-0.5 + t))*
              Sin(2*x) + 
              Sin(3*(-0.5 + t))*
              Sin(3*x)))*Sin(y) + 
         3*(7*Sin(1. - 2.*t)*Sin(x)*
             Sin(y) - 
            (7*Sin(4*(-0.5 + t))*
              (2*x*(1 + x)*
             Cos(2*x) + Sin(2*x)) + 
              2*Cos(2*(-0.5 + t))*
              (-7*x*(1 + x)*Cos(x)*
              Sin(0.5 - t) + 
              15*x*(1 + x)*Cos(3*x)*
              Sin(3*(-0.5 + t)) - 
              7*Sin(0.5 - t)*
             Sin(x) + 
              5*Sin(3*(-0.5 + t))*
              Sin(3*x)))*Sin(2*y)))\
       + 420*nu*
       (21*(1 + x)*Cos(x)*
          ((3 + x + 2*Power(x,2) + 
              Power(x,3))*
             Sin(1. - 2.*t) + 
            2*(6 + 5*x + 
              10*Power(x,2) + 
              5*Power(x,3))*
             Cos(2*(-0.5 + t))*
             Cos(y)*Sin(0.5 - t))*
          Sin(y) - 
         35*Cos(0.5 - t)*
          (2*(6 + 11*x + 
              15*Power(x,2) + 
              15*Power(x,3) + 
              5*Power(x,4))*Cos(2*x)*
             Sin(2*(-0.5 + t)) + 
            3*(3 + 8*x + 
              15*Power(x,2) + 
              15*Power(x,3) + 
              5*Power(x,4))*Cos(3*x)*
             Sin(3*(-0.5 + t)) + 
            7*Sin(2*(-0.5 + t))*
             Sin(2*x) + 
            26*x*Sin(2*(-0.5 + t))*
             Sin(2*x) + 
            13*Power(x,2)*
             Sin(2*(-0.5 + t))*
             Sin(2*x) + 
            11*Sin(3*(-0.5 + t))*
             Sin(3*x) + 
            28*x*Sin(3*(-0.5 + t))*
             Sin(3*x) + 
            14*Power(x,2)*
             Sin(3*(-0.5 + t))*
             Sin(3*x))*Sin(y) + 
         3*(7*(-1 + 4*x + 
              2*Power(x,2))*
             Sin(1. - 2.*t)*Sin(x)*
             Sin(y) - 
            (7*Sin(4*(-0.5 + t))*
              (2*
              (3 + 7*x + 
              12*Power(x,2) + 
              12*Power(x,3) + 
              4*Power(x,4))*Cos(2*x)\
              + (5 + 16*x + 
              8*Power(x,2))*Sin(2*x))
               + Cos(2*(-0.5 + t))*
              (15*
              (6 + 19*x + 
              39*Power(x,2) + 
              39*Power(x,3) + 
              13*Power(x,4))*
              Cos(3*x)*
              Sin(3*(-0.5 + t)) - 
              7*
              (1 + 14*x + 
              7*Power(x,2))*
              Sin(0.5 - t)*Sin(x) + 
              5*
              (25 + 62*x + 
              31*Power(x,2))*
              Sin(3*(-0.5 + t))*
              Sin(3*x)))*Sin(2*y)))\
       + (6*(-7*
             (5 + 5*x + 
              x*Cos(2*y)*
              Sin(4*(-0.5 + t))*
              Sin(2*x)) + 
            2*x*Cos(2*(-0.5 + t))*
             Cos(2*y)*
             (7*Sin(0.5 - t)*
             Sin(x) - 
              5*Sin(3*(-0.5 + t))*
              Sin(3*x))) + 
         7*x*Cos(0.5 - t)*Cos(y)*
          (6*Sin(0.5 - t)*Sin(x) - 
            5*(2*Sin(2*(-0.5 + t))*
              Sin(2*x) + 
              Sin(3*(-0.5 + t))*
              Sin(3*x))))*
       (-42*(1 + x)*Cos(x)*
          (Sin(1. - 2.*t) + 
            4*Cos(2*(-0.5 + t))*
             Cos(y)*Sin(0.5 - t))*
          Sin(y) - 
         35*Cos(0.5 - t)*
          (-8*(1 + x)*Cos(2*x)*
             Sin(2*(-0.5 + t)) - 
            6*(1 + x)*Cos(3*x)*
             Sin(3*(-0.5 + t)) + 
            4*Sin(2*(-0.5 + t))*
             Sin(2*x) + 
            8*x*Sin(2*(-0.5 + t))*
             Sin(2*x) + 
            16*Power(x,2)*
             Sin(2*(-0.5 + t))*
             Sin(2*x) + 
            8*Power(x,3)*
             Sin(2*(-0.5 + t))*
             Sin(2*x) + 
            2*Sin(3*(-0.5 + t))*
             Sin(3*x) + 
            9*x*Sin(3*(-0.5 + t))*
             Sin(3*x) + 
            18*Power(x,2)*
             Sin(3*(-0.5 + t))*
             Sin(3*x) + 
            9*Power(x,3)*
             Sin(3*(-0.5 + t))*
             Sin(3*x))*Sin(y) + 
         3*(7*(2 + x + 
              2*Power(x,2) + 
              Power(x,3))*
             Sin(1. - 2.*t)*Sin(x)*
             Sin(y) + 
            2*(-7*Sin(4*(-0.5 + t))*
              (-2*(1 + x)*Cos(2*x) + 
              (1 + 2*x + 
              4*Power(x,2) + 
              2*Power(x,3))*Sin(2*x))
               + Cos(2*(-0.5 + t))*
              (30*(1 + x)*Cos(3*x)*
              Sin(3*(-0.5 + t)) + 
              7*
              (2 + x + 
              2*Power(x,2) + 
              Power(x,3))*
              Sin(0.5 - t)*Sin(x) - 
              5*
              (2 + 9*x + 
              18*Power(x,2) + 
              9*Power(x,3))*
              Sin(3*(-0.5 + t))*
              Sin(3*x)))*Sin(2*y)))\
       - 210*x*Power(1 + x,3)*
       (42*Power(Cos(0.5 - t),2)*
          Cos(x)*Sin(y) + 
         140*Cos(2*x)*Sin(0.5 - t)*
          Sin(2*(-0.5 + t))*Sin(y) - 
         42*Cos(x)*Sin(0.5 - t)*
          (Sin(0.5 - t) - 
            4*Cos(y)*
             Sin(2*(-0.5 + t)))*
          Sin(y) + 
         105*Cos(3*x)*Sin(0.5 - t)*
          Sin(3*(-0.5 + t))*Sin(y) + 
         168*Power(Cos(2*(-0.5 + t)),
           2)*Cos(2*x)*Sin(2*y) + 
         270*Cos(2*(-0.5 + t))*
          Cos(3*(-0.5 + t))*Cos(3*x)*
          Sin(2*y) - 
         168*Cos(2*x)*
          Power(Sin(2*(-0.5 + t)),2)*
          Sin(2*y) - 
         180*Cos(3*x)*
          Sin(2*(-0.5 + t))*
          Sin(3*(-0.5 + t))*Sin(2*y)\
          + 7*Cos(0.5 - t)*
          (45*Cos(3*(-0.5 + t))*
             Cos(3*x)*Sin(y) + 
            Cos(2*(-0.5 + t))*
             (40*Cos(2*x)*Sin(y) + 
              6*Cos(x)*Sin(2*y)))) + 
      210*x*Power(1 + x,2)*
       (42*Power(Cos(0.5 - t),2)*
          Sin(x)*Sin(y) - 
         42*Power(Sin(0.5 - t),2)*
          Sin(x)*Sin(y) + 
         6*(14*
             Power(Cos(2*(-0.5 + t)),
              2)*Sin(2*x) + 
            15*Cos(2*(-0.5 + t))*
             Cos(3*(-0.5 + t))*
             Sin(3*x) - 
            2*Sin(2*(-0.5 + t))*
             (7*Sin(2*(-0.5 + t))*
              Sin(2*x) + 
              5*Sin(3*(-0.5 + t))*
              Sin(3*x)))*Sin(2*y) + 
         7*Sin(0.5 - t)*
          (5*Sin(3*(-0.5 + t))*
             Sin(3*x)*Sin(y) + 
            2*Sin(2*(-0.5 + t))*
             (5*Sin(2*x)*Sin(y) + 
              6*Sin(x)*Sin(2*y))) + 
         7*Cos(0.5 - t)*
          (15*Cos(3*(-0.5 + t))*
             Sin(3*x)*Sin(y) + 
            Cos(2*(-0.5 + t))*
             (20*Sin(2*x)*Sin(y) + 
              6*Sin(x)*Sin(2*y)))) - 
      210*Power(1 + x,3)*
       (42*Power(Cos(0.5 - t),2)*
          Sin(x)*Sin(y) - 
         42*Power(Sin(0.5 - t),2)*
          Sin(x)*Sin(y) + 
         6*(14*
             Power(Cos(2*(-0.5 + t)),
              2)*Sin(2*x) + 
            15*Cos(2*(-0.5 + t))*
             Cos(3*(-0.5 + t))*
             Sin(3*x) - 
            2*Sin(2*(-0.5 + t))*
             (7*Sin(2*(-0.5 + t))*
              Sin(2*x) + 
              5*Sin(3*(-0.5 + t))*
              Sin(3*x)))*Sin(2*y) + 
         7*Sin(0.5 - t)*
          (5*Sin(3*(-0.5 + t))*
             Sin(3*x)*Sin(y) + 
            2*Sin(2*(-0.5 + t))*
             (5*Sin(2*x)*Sin(y) + 
              6*Sin(x)*Sin(2*y))) + 
         7*Cos(0.5 - t)*
          (15*Cos(3*(-0.5 + t))*
             Sin(3*x)*Sin(y) + 
            Cos(2*(-0.5 + t))*
             (20*Sin(2*x)*Sin(y) + 
              6*Sin(x)*Sin(2*y)))))/
	      (44100.*Power(1 + x,4)));
}





int main(int argc, char** argv)
{
  try
    {
      Time timer("total");
      timer.start();
      int nx = 32;
      Sundance::setOption("nx", nx, "grid size");

      int nSteps = 32;
      Sundance::setOption("nSteps", nSteps, "number of timesteps");

      double tCur = 0.5; // time can't hit 0 + pi*k since q has a csc(t) term
      Sundance::setOption("iInit", tCur, "initial time");
      
      double tFinal = 2.0;
      Sundance::setOption("tFinal", tFinal, "final time");

      int order=2;
      Sundance::setOption("order", order, "stepper order");

      int verb = 1;
      Sundance::setOption("verb", verb, "verbosity");
      Sundance::init(&argc, &argv);

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedRectangleMesher. */
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedRectangleMesher(0.0, 2.0*M_PI, nx, 1,
                                                         0.0, M_PI, nx, 1,
                                                         meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter bdry = new BoundaryCellFilter();
      CellFilter top = bdry.subset( new TopPointTest() );
      CellFilter right = bdry.subset( new RightPointTest() );
      CellFilter bottom = bdry.subset( new BottomPointTest() );
      CellFilter left = bdry.subset( new LeftPointTest() );
      CellFilter nodes = new DimensionalCellFilter(0);
      CellFilter corner = nodes.subset(new CornerPointTest());

      /* Bases for velocity and pressure */
      BasisFamily vBas  = new Lagrange(2);
      BasisFamily pBas  = new Lagrange(1);

      /* Create the FlowState object */

      RCP<ProblemDescription> prob = rcp(new ProblemDescription(mesh));
      prob->addCondition(BodyForce, interior, StateDependentFunction(qExFunc));
      
      // Create No Penetration Boundary Conditions on top and bottom
      prob->addCondition(FreeSlip,top,StateDependentFunction(0.0));
      prob->addCondition(FreeSlip,bottom,StateDependentFunction(0.0));

      //Set inflow and outflow conditions
      prob->addCondition(NoSlip,left,StateDependentFunction(List(1.0,0.0)));
      prob->addCondition(Traction,right,StateDependentFunction(0.0));

      // Set IC
      prob->setIC(StateDependentFunction(List(1.0,0.0)),StateDependentFunction(pExFunc));
      prob->setNuEff(1.0);
      

      PCDControlParams pcdParams;
      NewtonControlParams newtonParams;

      SteppingType method = FullyImplicit;
      //SteppingType method = SemiImplicit;
      //pcdParams.FParams.method = "CG";

      PCDStepper stepper(prob, pcdParams, newtonParams, order, method);
      RCP<FlowState> state = stepper.state();

      double dt = (tFinal-tCur)/nSteps;
      
      state->initialize(tCur, dt, uExFunc, pExFunc);

      string outDir = "Results/";
   //   system( ("mkdir -p " + outDir).c_str() ); Done in DefaultOutputManager constructor
      string filename = outDir + "MMS_Inflow_nx" + Teuchos::toString(nx)
	+ "-nt-" + Teuchos::toString(nSteps);

      ErrorChecker check(state, uExFunc, pExFunc);

//Added by me to write out the snapshots
      RCP<DefaultOutputManager> output
	= rcp(new DefaultOutputManager(filename, "st", new ExodusWriterFactory()));

      // Same logic as VientoUniformStepController.cpp; I have taken it out so I can
      // utilize ErrorChecker
      for (int i=0; i<=nSteps; i++)
	{
	  Tabs tab1;
	  SUNDANCE_MSG1(verb, tab1<< "step " << i << " time=" << tCur);
	  if (i>0)
	    {
	      stepper.step(dt, verb);
	    }
	  Tabs tab2;
	  SUNDANCE_MSG1(verb, tab2<< "CFL=" << state->minHOverV() );
	  FieldWriter w = new VTKWriter(filename + "/" + Teuchos::toString(i));
	  check.write(w);
	  output->write(i, stepper.state(), 1);
	}

      timer.stop();
      Out::root() << "nx=" << nx << " nt=" << nSteps << " uErr=" << check.uErrL2()
		  << " pErr=" << check.pErrL2() << " runtime=" << timer.totalElapsedTime() << endl;

    }
  catch(std::exception& e)
    {
      Sundance::handleException(e);
    }
  Sundance::finalize(); return Sundance::testStatus(); 
}
