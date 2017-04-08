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

  return List(1 + x*Power(-2*Pi + x,2)*((Cos(t)*Cos(y)*Sin(t)*Sin(x))/5. + 
       (2*Cos(2*t)*Cos(2*y)*Sin(t)*Sin(x))/5. + 
       (Cos(t)*Cos(y)*Sin(2*t)*Sin(2*x))/3. + 
       (2*Cos(2*t)*Cos(2*y)*Sin(2*t)*Sin(2*x))/5. + 
       (Cos(t)*Cos(y)*Sin(3*t)*Sin(3*x))/6. + 
       (2*Cos(2*t)*Cos(2*y)*Sin(3*t)*Sin(3*x))/7.),
   -(x*Power(-2*Pi + x,2)*((Cos(t)*Cos(x)*Sin(t)*Sin(y))/5. + 
         (2*Cos(t)*Cos(2*x)*Sin(2*t)*Sin(y))/3. + 
         (Cos(t)*Cos(3*x)*Sin(3*t)*Sin(y))/2. + 
         (Cos(2*t)*Cos(x)*Sin(t)*Sin(2*y))/5. + 
         (2*Cos(2*t)*Cos(2*x)*Sin(2*t)*Sin(2*y))/5. + 
         (3*Cos(2*t)*Cos(3*x)*Sin(3*t)*Sin(2*y))/7.)) - 
    2*x*(-2*Pi + x)*((Cos(t)*Sin(t)*Sin(x)*Sin(y))/5. + 
       (Cos(t)*Sin(2*t)*Sin(2*x)*Sin(y))/3. + (Cos(t)*Sin(3*t)*Sin(3*x)*Sin(y))/6. + 
       (Cos(2*t)*Sin(t)*Sin(x)*Sin(2*y))/5. + 
       (Cos(2*t)*Sin(2*t)*Sin(2*x)*Sin(2*y))/5. + 
       (Cos(2*t)*Sin(3*t)*Sin(3*x)*Sin(2*y))/7.) - 
    Power(-2*Pi + x,2)*((Cos(t)*Sin(t)*Sin(x)*Sin(y))/5. + 
       (Cos(t)*Sin(2*t)*Sin(2*x)*Sin(y))/3. + (Cos(t)*Sin(3*t)*Sin(3*x)*Sin(y))/6. + 
       (Cos(2*t)*Sin(t)*Sin(x)*Sin(2*y))/5. + 
       (Cos(2*t)*Sin(2*t)*Sin(2*x)*Sin(2*y))/5. + 
			(Cos(2*t)*Sin(3*t)*Sin(3*x)*Sin(2*y))/7.));
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

  return List(-1 + Cos(x)*Cos(y) + x*Power(-2*Pi + x,2)*
     ((Power(Cos(t),2)*Cos(y)*Sin(x))/5. + (2*Cos(t)*Cos(2*t)*Cos(2*y)*Sin(x))/5. - 
       (Cos(y)*Power(Sin(t),2)*Sin(x))/5. - (4*Cos(2*y)*Sin(t)*Sin(2*t)*Sin(x))/5. + 
       (2*Cos(t)*Cos(2*t)*Cos(y)*Sin(2*x))/3. + 
       (4*Power(Cos(2*t),2)*Cos(2*y)*Sin(2*x))/5. - 
       (Cos(y)*Sin(t)*Sin(2*t)*Sin(2*x))/3. - 
       (4*Cos(2*y)*Power(Sin(2*t),2)*Sin(2*x))/5. + 
       (Cos(t)*Cos(3*t)*Cos(y)*Sin(3*x))/2. + 
       (6*Cos(2*t)*Cos(3*t)*Cos(2*y)*Sin(3*x))/7. - 
       (Cos(y)*Sin(t)*Sin(3*t)*Sin(3*x))/6. - 
       (4*Cos(2*y)*Sin(2*t)*Sin(3*t)*Sin(3*x))/7.) + 
    (x*Power(-2*Pi + x,2)*((Cos(t)*Cos(x)*Cos(y)*Sin(t))/5. + 
          (2*Cos(2*t)*Cos(x)*Cos(2*y)*Sin(t))/5. + 
          (2*Cos(t)*Cos(2*x)*Cos(y)*Sin(2*t))/3. + 
          (4*Cos(2*t)*Cos(2*x)*Cos(2*y)*Sin(2*t))/5. + 
          (Cos(t)*Cos(3*x)*Cos(y)*Sin(3*t))/2. + 
          (6*Cos(2*t)*Cos(3*x)*Cos(2*y)*Sin(3*t))/7.) + 
       2*x*(-2*Pi + x)*((Cos(t)*Cos(y)*Sin(t)*Sin(x))/5. + 
          (2*Cos(2*t)*Cos(2*y)*Sin(t)*Sin(x))/5. + 
          (Cos(t)*Cos(y)*Sin(2*t)*Sin(2*x))/3. + 
          (2*Cos(2*t)*Cos(2*y)*Sin(2*t)*Sin(2*x))/5. + 
          (Cos(t)*Cos(y)*Sin(3*t)*Sin(3*x))/6. + 
          (2*Cos(2*t)*Cos(2*y)*Sin(3*t)*Sin(3*x))/7.) + 
       Power(-2*Pi + x,2)*((Cos(t)*Cos(y)*Sin(t)*Sin(x))/5. + 
          (2*Cos(2*t)*Cos(2*y)*Sin(t)*Sin(x))/5. + 
          (Cos(t)*Cos(y)*Sin(2*t)*Sin(2*x))/3. + 
          (2*Cos(2*t)*Cos(2*y)*Sin(2*t)*Sin(2*x))/5. + 
          (Cos(t)*Cos(y)*Sin(3*t)*Sin(3*x))/6. + 
          (2*Cos(2*t)*Cos(2*y)*Sin(3*t)*Sin(3*x))/7.))*
     (1 + x*Power(-2*Pi + x,2)*((Cos(t)*Cos(y)*Sin(t)*Sin(x))/5. + 
          (2*Cos(2*t)*Cos(2*y)*Sin(t)*Sin(x))/5. + 
          (Cos(t)*Cos(y)*Sin(2*t)*Sin(2*x))/3. + 
          (2*Cos(2*t)*Cos(2*y)*Sin(2*t)*Sin(2*x))/5. + 
          (Cos(t)*Cos(y)*Sin(3*t)*Sin(3*x))/6. + 
          (2*Cos(2*t)*Cos(2*y)*Sin(3*t)*Sin(3*x))/7.)) - 
    nu*(2*(2*x*(-2*Pi + x) + Power(-2*Pi + x,2))*
        ((Cos(t)*Cos(x)*Cos(y)*Sin(t))/5. + (2*Cos(2*t)*Cos(x)*Cos(2*y)*Sin(t))/5. + 
          (2*Cos(t)*Cos(2*x)*Cos(y)*Sin(2*t))/3. + 
          (4*Cos(2*t)*Cos(2*x)*Cos(2*y)*Sin(2*t))/5. + 
          (Cos(t)*Cos(3*x)*Cos(y)*Sin(3*t))/2. + 
          (6*Cos(2*t)*Cos(3*x)*Cos(2*y)*Sin(3*t))/7.) + 
       x*Power(-2*Pi + x,2)*(-(Cos(t)*Cos(y)*Sin(t)*Sin(x))/5. - 
          (2*Cos(2*t)*Cos(2*y)*Sin(t)*Sin(x))/5. - 
          (4*Cos(t)*Cos(y)*Sin(2*t)*Sin(2*x))/3. - 
          (8*Cos(2*t)*Cos(2*y)*Sin(2*t)*Sin(2*x))/5. - 
          (3*Cos(t)*Cos(y)*Sin(3*t)*Sin(3*x))/2. - 
          (18*Cos(2*t)*Cos(2*y)*Sin(3*t)*Sin(3*x))/7.) + 
       x*Power(-2*Pi + x,2)*(-(Cos(t)*Cos(y)*Sin(t)*Sin(x))/5. - 
          (8*Cos(2*t)*Cos(2*y)*Sin(t)*Sin(x))/5. - 
          (Cos(t)*Cos(y)*Sin(2*t)*Sin(2*x))/3. - 
          (8*Cos(2*t)*Cos(2*y)*Sin(2*t)*Sin(2*x))/5. - 
          (Cos(t)*Cos(y)*Sin(3*t)*Sin(3*x))/6. - 
          (8*Cos(2*t)*Cos(2*y)*Sin(3*t)*Sin(3*x))/7.) + 
       (2*x + 4*(-2*Pi + x))*((Cos(t)*Cos(y)*Sin(t)*Sin(x))/5. + 
          (2*Cos(2*t)*Cos(2*y)*Sin(t)*Sin(x))/5. + 
          (Cos(t)*Cos(y)*Sin(2*t)*Sin(2*x))/3. + 
          (2*Cos(2*t)*Cos(2*y)*Sin(2*t)*Sin(2*x))/5. + 
          (Cos(t)*Cos(y)*Sin(3*t)*Sin(3*x))/6. + 
          (2*Cos(2*t)*Cos(2*y)*Sin(3*t)*Sin(3*x))/7.)) + 
    x*Power(-2*Pi + x,2)*(-(Cos(t)*Sin(t)*Sin(x)*Sin(y))/5. - 
       (Cos(t)*Sin(2*t)*Sin(2*x)*Sin(y))/3. - (Cos(t)*Sin(3*t)*Sin(3*x)*Sin(y))/6. - 
       (4*Cos(2*t)*Sin(t)*Sin(x)*Sin(2*y))/5. - 
       (4*Cos(2*t)*Sin(2*t)*Sin(2*x)*Sin(2*y))/5. - 
       (4*Cos(2*t)*Sin(3*t)*Sin(3*x)*Sin(2*y))/7.)*
     (-(x*Power(-2*Pi + x,2)*((Cos(t)*Cos(x)*Sin(t)*Sin(y))/5. + 
            (2*Cos(t)*Cos(2*x)*Sin(2*t)*Sin(y))/3. + 
            (Cos(t)*Cos(3*x)*Sin(3*t)*Sin(y))/2. + 
            (Cos(2*t)*Cos(x)*Sin(t)*Sin(2*y))/5. + 
            (2*Cos(2*t)*Cos(2*x)*Sin(2*t)*Sin(2*y))/5. + 
            (3*Cos(2*t)*Cos(3*x)*Sin(3*t)*Sin(2*y))/7.)) - 
       2*x*(-2*Pi + x)*((Cos(t)*Sin(t)*Sin(x)*Sin(y))/5. + 
          (Cos(t)*Sin(2*t)*Sin(2*x)*Sin(y))/3. + 
          (Cos(t)*Sin(3*t)*Sin(3*x)*Sin(y))/6. + 
          (Cos(2*t)*Sin(t)*Sin(x)*Sin(2*y))/5. + 
          (Cos(2*t)*Sin(2*t)*Sin(2*x)*Sin(2*y))/5. + 
          (Cos(2*t)*Sin(3*t)*Sin(3*x)*Sin(2*y))/7.) - 
       Power(-2*Pi + x,2)*((Cos(t)*Sin(t)*Sin(x)*Sin(y))/5. + 
          (Cos(t)*Sin(2*t)*Sin(2*x)*Sin(y))/3. + 
          (Cos(t)*Sin(3*t)*Sin(3*x)*Sin(y))/6. + 
          (Cos(2*t)*Sin(t)*Sin(x)*Sin(2*y))/5. + 
          (Cos(2*t)*Sin(2*t)*Sin(2*x)*Sin(2*y))/5. + 
          (Cos(2*t)*Sin(3*t)*Sin(3*x)*Sin(2*y))/7.)),
   -(Sin(x)*Sin(y)) - x*Power(-2*Pi + x,2)*
     ((Power(Cos(t),2)*Cos(x)*Sin(y))/5. + (4*Cos(t)*Cos(2*t)*Cos(2*x)*Sin(y))/3. + 
       (3*Cos(t)*Cos(3*t)*Cos(3*x)*Sin(y))/2. - (Cos(x)*Power(Sin(t),2)*Sin(y))/5. - 
       (2*Cos(2*x)*Sin(t)*Sin(2*t)*Sin(y))/3. - 
       (Cos(3*x)*Sin(t)*Sin(3*t)*Sin(y))/2. + (Cos(t)*Cos(2*t)*Cos(x)*Sin(2*y))/5. + 
       (4*Power(Cos(2*t),2)*Cos(2*x)*Sin(2*y))/5. + 
       (9*Cos(2*t)*Cos(3*t)*Cos(3*x)*Sin(2*y))/7. - 
       (2*Cos(x)*Sin(t)*Sin(2*t)*Sin(2*y))/5. - 
       (4*Cos(2*x)*Power(Sin(2*t),2)*Sin(2*y))/5. - 
       (6*Cos(3*x)*Sin(2*t)*Sin(3*t)*Sin(2*y))/7.) - 
    2*x*(-2*Pi + x)*((Power(Cos(t),2)*Sin(x)*Sin(y))/5. - 
       (Power(Sin(t),2)*Sin(x)*Sin(y))/5. + (2*Cos(t)*Cos(2*t)*Sin(2*x)*Sin(y))/3. - 
       (Sin(t)*Sin(2*t)*Sin(2*x)*Sin(y))/3. + (Cos(t)*Cos(3*t)*Sin(3*x)*Sin(y))/2. - 
       (Sin(t)*Sin(3*t)*Sin(3*x)*Sin(y))/6. + (Cos(t)*Cos(2*t)*Sin(x)*Sin(2*y))/5. - 
       (2*Sin(t)*Sin(2*t)*Sin(x)*Sin(2*y))/5. + 
       (2*Power(Cos(2*t),2)*Sin(2*x)*Sin(2*y))/5. - 
       (2*Power(Sin(2*t),2)*Sin(2*x)*Sin(2*y))/5. + 
       (3*Cos(2*t)*Cos(3*t)*Sin(3*x)*Sin(2*y))/7. - 
       (2*Sin(2*t)*Sin(3*t)*Sin(3*x)*Sin(2*y))/7.) - 
    Power(-2*Pi + x,2)*((Power(Cos(t),2)*Sin(x)*Sin(y))/5. - 
       (Power(Sin(t),2)*Sin(x)*Sin(y))/5. + (2*Cos(t)*Cos(2*t)*Sin(2*x)*Sin(y))/3. - 
       (Sin(t)*Sin(2*t)*Sin(2*x)*Sin(y))/3. + (Cos(t)*Cos(3*t)*Sin(3*x)*Sin(y))/2. - 
       (Sin(t)*Sin(3*t)*Sin(3*x)*Sin(y))/6. + (Cos(t)*Cos(2*t)*Sin(x)*Sin(2*y))/5. - 
       (2*Sin(t)*Sin(2*t)*Sin(x)*Sin(2*y))/5. + 
       (2*Power(Cos(2*t),2)*Sin(2*x)*Sin(2*y))/5. - 
       (2*Power(Sin(2*t),2)*Sin(2*x)*Sin(2*y))/5. + 
       (3*Cos(2*t)*Cos(3*t)*Sin(3*x)*Sin(2*y))/7. - 
       (2*Sin(2*t)*Sin(3*t)*Sin(3*x)*Sin(2*y))/7.) + 
    (1 + x*Power(-2*Pi + x,2)*((Cos(t)*Cos(y)*Sin(t)*Sin(x))/5. + 
          (2*Cos(2*t)*Cos(2*y)*Sin(t)*Sin(x))/5. + 
          (Cos(t)*Cos(y)*Sin(2*t)*Sin(2*x))/3. + 
          (2*Cos(2*t)*Cos(2*y)*Sin(2*t)*Sin(2*x))/5. + 
          (Cos(t)*Cos(y)*Sin(3*t)*Sin(3*x))/6. + 
          (2*Cos(2*t)*Cos(2*y)*Sin(3*t)*Sin(3*x))/7.))*
     (-4*x*(-2*Pi + x)*((Cos(t)*Cos(x)*Sin(t)*Sin(y))/5. + 
          (2*Cos(t)*Cos(2*x)*Sin(2*t)*Sin(y))/3. + 
          (Cos(t)*Cos(3*x)*Sin(3*t)*Sin(y))/2. + 
          (Cos(2*t)*Cos(x)*Sin(t)*Sin(2*y))/5. + 
          (2*Cos(2*t)*Cos(2*x)*Sin(2*t)*Sin(2*y))/5. + 
          (3*Cos(2*t)*Cos(3*x)*Sin(3*t)*Sin(2*y))/7.) - 
       2*Power(-2*Pi + x,2)*((Cos(t)*Cos(x)*Sin(t)*Sin(y))/5. + 
          (2*Cos(t)*Cos(2*x)*Sin(2*t)*Sin(y))/3. + 
          (Cos(t)*Cos(3*x)*Sin(3*t)*Sin(y))/2. + 
          (Cos(2*t)*Cos(x)*Sin(t)*Sin(2*y))/5. + 
          (2*Cos(2*t)*Cos(2*x)*Sin(2*t)*Sin(2*y))/5. + 
          (3*Cos(2*t)*Cos(3*x)*Sin(3*t)*Sin(2*y))/7.) - 
       x*Power(-2*Pi + x,2)*(-(Cos(t)*Sin(t)*Sin(x)*Sin(y))/5. - 
          (4*Cos(t)*Sin(2*t)*Sin(2*x)*Sin(y))/3. - 
          (3*Cos(t)*Sin(3*t)*Sin(3*x)*Sin(y))/2. - 
          (Cos(2*t)*Sin(t)*Sin(x)*Sin(2*y))/5. - 
          (4*Cos(2*t)*Sin(2*t)*Sin(2*x)*Sin(2*y))/5. - 
          (9*Cos(2*t)*Sin(3*t)*Sin(3*x)*Sin(2*y))/7.) - 
       2*x*((Cos(t)*Sin(t)*Sin(x)*Sin(y))/5. + 
          (Cos(t)*Sin(2*t)*Sin(2*x)*Sin(y))/3. + 
          (Cos(t)*Sin(3*t)*Sin(3*x)*Sin(y))/6. + 
          (Cos(2*t)*Sin(t)*Sin(x)*Sin(2*y))/5. + 
          (Cos(2*t)*Sin(2*t)*Sin(2*x)*Sin(2*y))/5. + 
          (Cos(2*t)*Sin(3*t)*Sin(3*x)*Sin(2*y))/7.) - 
       4*(-2*Pi + x)*((Cos(t)*Sin(t)*Sin(x)*Sin(y))/5. + 
          (Cos(t)*Sin(2*t)*Sin(2*x)*Sin(y))/3. + 
          (Cos(t)*Sin(3*t)*Sin(3*x)*Sin(y))/6. + 
          (Cos(2*t)*Sin(t)*Sin(x)*Sin(2*y))/5. + 
          (Cos(2*t)*Sin(2*t)*Sin(2*x)*Sin(2*y))/5. + 
          (Cos(2*t)*Sin(3*t)*Sin(3*x)*Sin(2*y))/7.)) + 
    (-(x*Power(-2*Pi + x,2)*((Cos(t)*Cos(x)*Cos(y)*Sin(t))/5. + 
            (2*Cos(2*t)*Cos(x)*Cos(2*y)*Sin(t))/5. + 
            (2*Cos(t)*Cos(2*x)*Cos(y)*Sin(2*t))/3. + 
            (4*Cos(2*t)*Cos(2*x)*Cos(2*y)*Sin(2*t))/5. + 
            (Cos(t)*Cos(3*x)*Cos(y)*Sin(3*t))/2. + 
            (6*Cos(2*t)*Cos(3*x)*Cos(2*y)*Sin(3*t))/7.)) - 
       2*x*(-2*Pi + x)*((Cos(t)*Cos(y)*Sin(t)*Sin(x))/5. + 
          (2*Cos(2*t)*Cos(2*y)*Sin(t)*Sin(x))/5. + 
          (Cos(t)*Cos(y)*Sin(2*t)*Sin(2*x))/3. + 
          (2*Cos(2*t)*Cos(2*y)*Sin(2*t)*Sin(2*x))/5. + 
          (Cos(t)*Cos(y)*Sin(3*t)*Sin(3*x))/6. + 
          (2*Cos(2*t)*Cos(2*y)*Sin(3*t)*Sin(3*x))/7.) - 
       Power(-2*Pi + x,2)*((Cos(t)*Cos(y)*Sin(t)*Sin(x))/5. + 
          (2*Cos(2*t)*Cos(2*y)*Sin(t)*Sin(x))/5. + 
          (Cos(t)*Cos(y)*Sin(2*t)*Sin(2*x))/3. + 
          (2*Cos(2*t)*Cos(2*y)*Sin(2*t)*Sin(2*x))/5. + 
          (Cos(t)*Cos(y)*Sin(3*t)*Sin(3*x))/6. + 
          (2*Cos(2*t)*Cos(2*y)*Sin(3*t)*Sin(3*x))/7.))*
     (-(x*Power(-2*Pi + x,2)*((Cos(t)*Cos(x)*Sin(t)*Sin(y))/5. + 
            (2*Cos(t)*Cos(2*x)*Sin(2*t)*Sin(y))/3. + 
            (Cos(t)*Cos(3*x)*Sin(3*t)*Sin(y))/2. + 
            (Cos(2*t)*Cos(x)*Sin(t)*Sin(2*y))/5. + 
            (2*Cos(2*t)*Cos(2*x)*Sin(2*t)*Sin(2*y))/5. + 
            (3*Cos(2*t)*Cos(3*x)*Sin(3*t)*Sin(2*y))/7.)) - 
       2*x*(-2*Pi + x)*((Cos(t)*Sin(t)*Sin(x)*Sin(y))/5. + 
          (Cos(t)*Sin(2*t)*Sin(2*x)*Sin(y))/3. + 
          (Cos(t)*Sin(3*t)*Sin(3*x)*Sin(y))/6. + 
          (Cos(2*t)*Sin(t)*Sin(x)*Sin(2*y))/5. + 
          (Cos(2*t)*Sin(2*t)*Sin(2*x)*Sin(2*y))/5. + 
          (Cos(2*t)*Sin(3*t)*Sin(3*x)*Sin(2*y))/7.) - 
       Power(-2*Pi + x,2)*((Cos(t)*Sin(t)*Sin(x)*Sin(y))/5. + 
          (Cos(t)*Sin(2*t)*Sin(2*x)*Sin(y))/3. + 
          (Cos(t)*Sin(3*t)*Sin(3*x)*Sin(y))/6. + 
          (Cos(2*t)*Sin(t)*Sin(x)*Sin(2*y))/5. + 
          (Cos(2*t)*Sin(2*t)*Sin(2*x)*Sin(2*y))/5. + 
          (Cos(2*t)*Sin(3*t)*Sin(3*x)*Sin(2*y))/7.)) - 
    nu*(-(x*Power(-2*Pi + x,2)*(-(Cos(t)*Cos(x)*Sin(t)*Sin(y))/5. - 
            (2*Cos(t)*Cos(2*x)*Sin(2*t)*Sin(y))/3. - 
            (Cos(t)*Cos(3*x)*Sin(3*t)*Sin(y))/2. - 
            (4*Cos(2*t)*Cos(x)*Sin(t)*Sin(2*y))/5. - 
            (8*Cos(2*t)*Cos(2*x)*Sin(2*t)*Sin(2*y))/5. - 
            (12*Cos(2*t)*Cos(3*x)*Sin(3*t)*Sin(2*y))/7.)) - 
       4*(-2*Pi + x)*((Cos(t)*Cos(x)*Sin(t)*Sin(y))/5. + 
          (2*Cos(t)*Cos(2*x)*Sin(2*t)*Sin(y))/3. + 
          (Cos(t)*Cos(3*x)*Sin(3*t)*Sin(y))/2. + 
          (Cos(2*t)*Cos(x)*Sin(t)*Sin(2*y))/5. + 
          (2*Cos(2*t)*Cos(2*x)*Sin(2*t)*Sin(2*y))/5. + 
          (3*Cos(2*t)*Cos(3*x)*Sin(3*t)*Sin(2*y))/7.) - 
       Power(-2*Pi + x,2)*(-(Cos(t)*Sin(t)*Sin(x)*Sin(y))/5. - 
          (4*Cos(t)*Sin(2*t)*Sin(2*x)*Sin(y))/3. - 
          (3*Cos(t)*Sin(3*t)*Sin(3*x)*Sin(y))/2. - 
          (Cos(2*t)*Sin(t)*Sin(x)*Sin(2*y))/5. - 
          (4*Cos(2*t)*Sin(2*t)*Sin(2*x)*Sin(2*y))/5. - 
          (9*Cos(2*t)*Sin(3*t)*Sin(3*x)*Sin(2*y))/7.) - 
       2*x*(-2*Pi + x)*(-(Cos(t)*Sin(t)*Sin(x)*Sin(y))/5. - 
          (Cos(t)*Sin(2*t)*Sin(2*x)*Sin(y))/3. - 
          (Cos(t)*Sin(3*t)*Sin(3*x)*Sin(y))/6. - 
          (4*Cos(2*t)*Sin(t)*Sin(x)*Sin(2*y))/5. - 
          (4*Cos(2*t)*Sin(2*t)*Sin(2*x)*Sin(2*y))/5. - 
          (4*Cos(2*t)*Sin(3*t)*Sin(3*x)*Sin(2*y))/7.) - 
       Power(-2*Pi + x,2)*(-(Cos(t)*Sin(t)*Sin(x)*Sin(y))/5. - 
          (Cos(t)*Sin(2*t)*Sin(2*x)*Sin(y))/3. - 
          (Cos(t)*Sin(3*t)*Sin(3*x)*Sin(y))/6. - 
          (4*Cos(2*t)*Sin(t)*Sin(x)*Sin(2*y))/5. - 
          (4*Cos(2*t)*Sin(2*t)*Sin(2*x)*Sin(2*y))/5. - 
          (4*Cos(2*t)*Sin(3*t)*Sin(3*x)*Sin(2*y))/7.) - 
       2*((Cos(t)*Sin(t)*Sin(x)*Sin(y))/5. + (Cos(t)*Sin(2*t)*Sin(2*x)*Sin(y))/3. + 
          (Cos(t)*Sin(3*t)*Sin(3*x)*Sin(y))/6. + 
          (Cos(2*t)*Sin(t)*Sin(x)*Sin(2*y))/5. + 
          (Cos(2*t)*Sin(2*t)*Sin(2*x)*Sin(2*y))/5. + 
          (Cos(2*t)*Sin(3*t)*Sin(3*x)*Sin(2*y))/7.) - 
       4*((Cos(t)*Sin(t)*Sin(x)*Sin(y))/5. + (Cos(t)*Sin(2*t)*Sin(2*x)*Sin(y))/3. + 
          (Cos(t)*Sin(3*t)*Sin(3*x)*Sin(y))/6. + 
          (Cos(2*t)*Sin(t)*Sin(x)*Sin(2*y))/5. + 
          (Cos(2*t)*Sin(2*t)*Sin(2*x)*Sin(2*y))/5. + 
          (Cos(2*t)*Sin(3*t)*Sin(3*x)*Sin(2*y))/7. + 
          (-2*Pi + x)*((Cos(t)*Cos(x)*Sin(t)*Sin(y))/5. + 
             (2*Cos(t)*Cos(2*x)*Sin(2*t)*Sin(y))/3. + 
             (Cos(t)*Cos(3*x)*Sin(3*t)*Sin(y))/2. + 
             (Cos(2*t)*Cos(x)*Sin(t)*Sin(2*y))/5. + 
             (2*Cos(2*t)*Cos(2*x)*Sin(2*t)*Sin(2*y))/5. + 
             (3*Cos(2*t)*Cos(3*x)*Sin(3*t)*Sin(2*y))/7.)) - 
       2*x*(2*((Cos(t)*Cos(x)*Sin(t)*Sin(y))/5. + 
             (2*Cos(t)*Cos(2*x)*Sin(2*t)*Sin(y))/3. + 
             (Cos(t)*Cos(3*x)*Sin(3*t)*Sin(y))/2. + 
             (Cos(2*t)*Cos(x)*Sin(t)*Sin(2*y))/5. + 
             (2*Cos(2*t)*Cos(2*x)*Sin(2*t)*Sin(2*y))/5. + 
             (3*Cos(2*t)*Cos(3*x)*Sin(3*t)*Sin(2*y))/7.) + 
          (-2*Pi + x)*(-(Cos(t)*Sin(t)*Sin(x)*Sin(y))/5. - 
             (4*Cos(t)*Sin(2*t)*Sin(2*x)*Sin(y))/3. - 
             (3*Cos(t)*Sin(3*t)*Sin(3*x)*Sin(y))/2. - 
             (Cos(2*t)*Sin(t)*Sin(x)*Sin(2*y))/5. - 
             (4*Cos(2*t)*Sin(2*t)*Sin(2*x)*Sin(2*y))/5. - 
             (9*Cos(2*t)*Sin(3*t)*Sin(3*x)*Sin(2*y))/7.)) - 
       x*(Power(-2*Pi + x,2)*(-(Cos(t)*Cos(x)*Sin(t)*Sin(y))/5. - 
             (8*Cos(t)*Cos(2*x)*Sin(2*t)*Sin(y))/3. - 
             (9*Cos(t)*Cos(3*x)*Sin(3*t)*Sin(y))/2. - 
             (Cos(2*t)*Cos(x)*Sin(t)*Sin(2*y))/5. - 
             (8*Cos(2*t)*Cos(2*x)*Sin(2*t)*Sin(2*y))/5. - 
             (27*Cos(2*t)*Cos(3*x)*Sin(3*t)*Sin(2*y))/7.) + 
          2*((Cos(t)*Cos(x)*Sin(t)*Sin(y))/5. + 
             (2*Cos(t)*Cos(2*x)*Sin(2*t)*Sin(y))/3. + 
             (Cos(t)*Cos(3*x)*Sin(3*t)*Sin(y))/2. + 
             (Cos(2*t)*Cos(x)*Sin(t)*Sin(2*y))/5. + 
             (2*Cos(2*t)*Cos(2*x)*Sin(2*t)*Sin(2*y))/5. + 
             (3*Cos(2*t)*Cos(3*x)*Sin(3*t)*Sin(2*y))/7.) + 
          4*(-2*Pi + x)*(-(Cos(t)*Sin(t)*Sin(x)*Sin(y))/5. - 
             (4*Cos(t)*Sin(2*t)*Sin(2*x)*Sin(y))/3. - 
             (3*Cos(t)*Sin(3*t)*Sin(3*x)*Sin(y))/2. - 
             (Cos(2*t)*Sin(t)*Sin(x)*Sin(2*y))/5. - 
             (4*Cos(2*t)*Sin(2*t)*Sin(2*x)*Sin(2*y))/5. - 
             (9*Cos(2*t)*Sin(3*t)*Sin(3*x)*Sin(2*y))/7.)) + 
       2*(-2*(-2*Pi + x)*((Cos(t)*Cos(x)*Sin(t)*Sin(y))/5. + 
             (2*Cos(t)*Cos(2*x)*Sin(2*t)*Sin(y))/3. + 
             (Cos(t)*Cos(3*x)*Sin(3*t)*Sin(y))/2. + 
             (Cos(2*t)*Cos(x)*Sin(t)*Sin(2*y))/5. + 
             (2*Cos(2*t)*Cos(2*x)*Sin(2*t)*Sin(2*y))/5. + 
             (3*Cos(2*t)*Cos(3*x)*Sin(3*t)*Sin(2*y))/7.) - 
          Power(-2*Pi + x,2)*(-(Cos(t)*Sin(t)*Sin(x)*Sin(y))/5. - 
             (4*Cos(t)*Sin(2*t)*Sin(2*x)*Sin(y))/3. - 
             (3*Cos(t)*Sin(3*t)*Sin(3*x)*Sin(y))/2. - 
             (Cos(2*t)*Sin(t)*Sin(x)*Sin(2*y))/5. - 
             (4*Cos(2*t)*Sin(2*t)*Sin(2*x)*Sin(2*y))/5. - 
			      (9*Cos(2*t)*Sin(3*t)*Sin(3*x)*Sin(2*y))/7.))));
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

      double tCur = 0.0; // time can't hit 0 + pi*k since q has a csc(t) term
      Sundance::setOption("iInit", tCur, "initial time");
      
      double tFinal = 2.0;
      Sundance::setOption("tFinal", tFinal, "final time");

      int order=2;
      Sundance::setOption("order", order, "stepper order");

      int verb = 1;
      Sundance::setOption("verbosity", verb, "verbosity");
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
		  << " pErr=" << check.pErrL2() << " runtime=" << timer.totalElapsedTime()
		  << endl << endl << endl;

    }
  catch(std::exception& e)
    {
      Sundance::handleException(e);
    }
  Sundance::finalize(); return Sundance::testStatus(); 
}
