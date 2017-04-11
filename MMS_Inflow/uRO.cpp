#include "Teuchos_GlobalMPISession.hpp"
#include "PlayaSerialVectorType.hpp"
#include "PlayaDenseSerialMatrix.hpp"
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaSerialEpetraAdapter.hpp" // Allows me to convert between Epetra and Serial Vectors
#include "PlayaOut.hpp"
#include <vector>


#include "VientoSnapshotIO.hpp" //For snapshotToMatrix
#include "denseSerialMatrixIO.hpp"
#include "QuadraticODERHSBase.hpp"
#include "MathematicaConverter.hpp"
#include "MMSQuadODE.hpp"
#include "MyNLO.hpp"
#include "velocityROM.hpp"

//For using newton-armijo
#include "PlayaDenseLUSolver.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "PlayaNewtonArmijoSolverImpl.hpp"
#include "PlayaNonlinearSolver.hpp"
#include "PlayaNonlinearSolverBuilder.hpp"


//Local files
#include "PlayaSVD.hpp"

#include "Sundance.hpp"

using std::vector;
using std::cout;
//using std::endl;
using namespace Teuchos;
using namespace Playa;
using namespace PlayaExprTemplates;




int main(int argc, char *argv[]) 
{
  try
    {
      Time timer("total");
      timer.start();
      
      int verbosity = 1;	
      Sundance::setOption("verbosity", verbosity, "verbosity level");
      
      int nx = 32;
      Sundance::setOption("nx", nx, "Number of elements along each axis");

      int nSteps = 32;
      Sundance::setOption("nSteps", nSteps, "Number of time steps");

      double tInit = 0.0; 
      Sundance::setOption("tInit", tInit, "initial time");

      double tFinal = 2.0;
      Sundance::setOption("tFinal", tFinal, "final time");

      bool AreMatrixAndTensorInFile = false;
      Sundance::setOption("MatrixAndTensorInFile", "MatrixAndTensorNotInFile", AreMatrixAndTensorInFile, "true if the matrix and tensor are available text files. false if they need to be created");

      Sundance::init(&argc, &argv);
      GlobalMPISession session(&argc, &argv);

      // Define our coordinates
      Expr x = new CoordExpr(0,"x");
      Expr y = new CoordExpr(1,"y");
      Expr t = new Sundance::Parameter(tInit);
      double deltat = (tFinal-tInit)/nSteps;
      double nu = 1.0;

      // Define initial conditions, i.e. u(x,y,0)
      Expr u0 = List(1.0,0.0);

      // Define the forcing term
      Expr q = List(-1 + Cos(x)*Cos(y) + x*Power(-2*Pi + x,2)*
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


      // Define our solution
      Expr uExact = List(1 + x*Power(-2*Pi + x,2)*((Cos(t)*Cos(y)*Sin(t)*Sin(x))/5. + 
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

      Expr pExact = 2*Pi - x + Cos(y)*Sin(x);


      

      // Define our mesh
      MeshType meshType = new Sundance::BasicSimplicialMeshType();
      double xmin = 0.0;
      double xmax = 2.0*Pi;
      double ymin = 0.0;
      double ymax = Pi;
      MeshSource mesher = new Sundance::PartitionedRectangleMesher(xmin,xmax,nx,1,ymin,ymax,nx,1,meshType,0);
      Mesh mesh = mesher.getMesh();



      // Read the snapshots into a matrix
      //      string outDir = "/home/sirush/PhDResearch/ODETest/Results";
      string outDir = "Results/MMS_Inflow_nx" + Teuchos::toString(nx)
	+ "-nt-" + Teuchos::toString(nSteps);
      string tag = "st-v";
      string filename = outDir + "/" + tag;

      // Create a BasisFamily to express our solution
      Array<Sundance::BasisFamily> ubasis = List(new Sundance::Lagrange(2), new Sundance::Lagrange(2)); // 2nd order Piece-Wise Quad Lagrange in 2D
      // Define our vector type
      Playa::VectorType<double> epetraVecType = new Playa::EpetraVectorType();
      Sundance::DiscreteSpace velocityDS(mesh, ubasis, epetraVecType);

      string NLParamFile = "playa-newton-armijo.xml";

      // Create a velocityROM object
      velocityROM ROM(filename, NLParamFile, velocityDS, u0, q, t, nSteps, deltat,
		      .99999, verbosity);
      ROM.initialize();
      ROM.generate_alpha();
      Array<Expr> uRO(ROM.get_uRO() );
     
      VectorType<double> time_vecType = new SerialVectorType();
      VectorSpace<double> time_vecSpace = time_vecType.createEvenlyPartitionedSpace(MPIComm::self(), nSteps+1.0);

      //Needed for the integral
      CellFilter interior = new MaximalCellFilter();
      QuadratureFamily quad = new GaussianQuadrature(6);


      // Compute alpha "exactly"
      Array<Expr> phi(ROM.get_phi() ); // These are the POD basis functions
      



      // Based off the value for R, create an appropriate VectorSpace<double>
      int R = phi.length();
      VectorType<double> R_vecType = new SerialVectorType();
      VectorSpace<double> R_vecSpace = R_vecType.createEvenlyPartitionedSpace(MPIComm::self(), R);
      
      //Find the exact alphas
      Array<Vector<double> > alpha(nSteps+1);
      for(int count = 0; count<alpha.length(); count++)
	alpha[count] = R_vecSpace.createMember();



      // Calculate ubar(x)
      LinearOperator<double> Wprime = snapshotToMatrix(filename, nSteps, mesh);
      SUNDANCE_ROOT_MSG2(verbosity, "Size of W: " << Wprime.range().dim() << " by " << Wprime.domain().dim());
      Vector<double> ubarVec = Wprime.range().createMember();
      Vector<double> ones = Wprime.domain().createMember();
      ones.setToConstant(1.0);
      Wprime.apply(ones, ubarVec);
      ubarVec *= (1.0/ (nSteps+1.0) );
      Expr ubar = new DiscreteFunction(velocityDS, serialToEpetra(ubarVec));

      //Find the exact alphas
      for(int tIndex=0; tIndex < alpha.length(); tIndex++)
	{
	  t.setParameterValue(tInit+tIndex*deltat);
	  for(int r=0; r<R; r++)
	    {
	      // alpha_r(t_m) = ( uEx(t_m, x, y), phi[r] )
	      //FunctionalEvaluator ExactEvaluator(mesh, Integral(interior, (uExact-ubar)*phi[r], quad));
	      FunctionalEvaluator ExactEvaluator = FunctionalEvaluator(mesh, Integral(interior, (uExact-ubar)*phi[r], quad));
	      alpha[tIndex][r] = ExactEvaluator.evaluate();
	    }
	}

      Array<Vector<double> > soln(ROM.get_alpha() ); // Approximate alphas


      
      Vector<double> alphaError = time_vecSpace.createMember();
      for(int i=0; i < alpha.length(); i++)
	{
	  alphaError[i] = (alpha[i] - soln[i]).norm2();
	  if(verbosity>=2)
	    cout << "Error for alpha(t=" << i << "): " << alphaError[i] << endl;

	  if(verbosity>=3)
	    {
	      cout << "Exact alpha(t=" << i << "): " << endl << alpha[i] << endl;	
	      cout << "Approximate alpha(t=" << i << "): " << endl << soln[i] << endl << endl;
	    }
	}

      cout << "Run for nx = " << nx << ", nSteps = " << nSteps << endl;
      cout << "||alphaExact - alphaApprox||_2  :\t "  << alphaError.norm2() << endl;
      cout << "||alphaExact - alphaApprox||_inf:\t " << alphaError.normInf() << endl;
      

      SUNDANCE_ROOT_MSG2(verbosity, "Comparing uExact(t_n) to uRO(t_n)");
      Vector<double> l2norm = time_vecSpace.createMember();
      
      for(int time = 0; time < uRO.length(); time++)
	{
	  t.setParameterValue(tInit+time*deltat);
	  l2norm[time] = L2Norm(mesh, interior, uExact - uRO[time], quad);
	  double uNorm = L2Norm(mesh, interior, uExact, quad);
	  /* print the relative error */
	  SUNDANCE_ROOT_MSG2(verbosity, "Error for uRO at time " << time*deltat
			     << " = " << l2norm[time]/uNorm);
	}
      
      timer.stop();
  SUNDANCE_ROOT_MSG1(verbosity, "Number of velocity modes kept: " + Teuchos::toString(ROM.get_phi().size()));
      Out::root() << "||uExact - uRO||_2  :\t " << l2norm.norm2() << endl;
      Out::root() << "||uExact - uRO||_inf:\t " << l2norm.normInf() << endl;
      Out::root() << "runtime=" << timer.totalElapsedTime() << endl << endl;
     
      // Visualize the results
      SUNDANCE_ROOT_MSG1(verbosity, "Writing results to file");
      string vtkDir = "Results/Visuals/uRO/";
      string vtkfilename = "nx"+Teuchos::toString(nx)+"nt"+Teuchos::toString(nSteps);
      vtkDir = vtkDir + vtkfilename + "/";
      system( ("mkdir -p " + vtkDir).c_str() ); 

      DiscreteSpace scalarDS(mesh, new Lagrange(1), new EpetraVectorType());
      L2Projector projector(velocityDS, uExact);
      //Write uExact for all the time steps
      for(int time=0; time< nSteps+1; time++)
	{
	  FieldWriter writer = new VTKWriter(vtkDir+vtkfilename+"step"+Teuchos::toString(time));
	  writer.addMesh(mesh);
	  t.setParameterValue(time*deltat);
	  L2Projector projectorRO(velocityDS, uRO[time]);
	  L2Projector uErrorProjector(velocityDS, uExact - uRO[time]);
	  Expr absErr = sqrt( (uExact - uRO[time])*(uExact - uRO[time]));
	  Expr absU = sqrt(uExact * uExact);
	  L2Projector uMagProj(scalarDS, absU);
	  L2Projector absErrorProj(scalarDS, absErr);
	  L2Projector relErrorProj(scalarDS, absErr / (absU + 1.0));
	  writer.addField("uMag", new ExprFieldWrapper(uMagProj.project()[0]) );
	  writer.addField("errAbs", new ExprFieldWrapper(absErrorProj.project()[0]) );
	  writer.addField("errRel", new ExprFieldWrapper(relErrorProj.project()[0]) );
	  writer.addField("uExact[0]", new ExprFieldWrapper(projector.project()[0]) );
	  writer.addField("uExact[1]", new ExprFieldWrapper(projector.project()[1]) );
	  writer.addField("uRO[0]", new ExprFieldWrapper(projectorRO.project()[0]) );
	  writer.addField("uRO[1]", new ExprFieldWrapper(projectorRO.project()[1]) );
	  writer.addField("uError[0]", new ExprFieldWrapper(uErrorProjector.project()[0]) );
	  writer.addField("uError[1]", new ExprFieldWrapper(uErrorProjector.project()[1]) );
      	  writer.write();	  
	}

      

	
    }
  catch(std::exception& e)
    {
      Out::root() << "Caught exception: " << e.what() << std::endl;
      return -1;
    }
}

