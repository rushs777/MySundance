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
#include "NSEProjectedODE.hpp"
#include "MyNLO.hpp"
#include "velocityROM.hpp"

//For using newton-armijo
#include "PlayaDenseLUSolver.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "PlayaNewtonArmijoSolverImpl.hpp"
#include "PlayaNonlinearSolver.hpp"
#include "PlayaNonlinearSolverBuilder.hpp"

#include "Sundance.hpp"

// Needed to write only significant digits of tol as part of the file name
#include <sstream>
#include <iomanip>



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
      
      int nx = 25;
      Sundance::setOption("nx", nx, "Number of elements along each axis");

      double tol = .999;
      Sundance::setOption("tol", tol, "Tolerance requirement for the number of basis functions to keep");

      int precision = 3;
      Sundance::setOption("precision",precision,"Number of significant digits to keep in a filename using tol");  
      
      int nSteps = 25;
      Sundance::setOption("nSteps", nSteps, "Number of time steps");

      double tInit = 0.0; 
      Sundance::setOption("tInit", tInit, "initial time");

      double tFinal = 1.0;
      Sundance::setOption("tFinal", tFinal, "final time");

      double Re = 1.0;
      Sundance::setOption("Re", Re, "Reynolds number");

      double U0 = 1.0;
      Sundance::setOption("U0", U0, "Characteristic velocity used for calculating nu");
  
      double L0 = 1.0;
      Sundance::setOption("L0", L0, "Characteristic length used for calculating nu");       

      bool AreMatrixAndTensorInFile = false;
      Sundance::setOption("MatrixAndTensorInFile", "MatrixAndTensorNotInFile", AreMatrixAndTensorInFile, "true if the matrix and tensor are available text files. false if they need to be created");

      Sundance::init(&argc, &argv);
      GlobalMPISession session(&argc, &argv);

      // Define our coordinates
      Expr x = new CoordExpr(0,"x");
      Expr y = new CoordExpr(1,"y");
      Expr t = new Sundance::Parameter(tInit);
      double deltat = (tFinal-tInit)/nSteps;
      double R = Re;

      // Define initial conditions, i.e. u(x,y,0)
      Expr u0 = List(1.0,0.0);

      // Define the forcing term
      Expr q = List((Pi*(-((-1 + Power(x,2))*(32*Pi*Cos(4*Pi*y)*Sin((Pi*t)/2.)*Power(Sin(Pi*x),2) +270*Pi*Cos(6*Pi*y)*Sin(Pi*t)*Power(Sin(2*Pi*x),2))*(4*Sin((Pi*t)/2.)*Sin(Pi*x)*(Pi*(-1 + Power(x,2))*Cos(Pi*x) +x*Sin(Pi*x))*Power(Sin(2*Pi*y),2)+ 15*Sin(Pi*t)*Sin(2*Pi*x)*(2*Pi*(-1 + Power(x,2))*Cos(2*Pi*x) +x*Sin(2*Pi*x))*Power(Sin(3*Pi*y),2))) -100*Pi*(-1 + Power(x,2))*(4*Cos((Pi*t)/2.)*Power(Sin(Pi*x),2)*Sin(4*Pi*y) +45*Cos(Pi*t)*Power(Sin(2*Pi*x),2)*Sin(6*Pi*y)) +(200*(8*Power(Pi,2)*(-1 + Power(x,2))*Power(Cos(Pi*x),2)*Sin((Pi*t)/2.)*Sin(4*Pi*y) -8*Sin((Pi*t)/2.)*Sin(Pi*x)*(-4*Pi*x*Cos(Pi*x) +(-1 +9*Power(Pi,2)*(-1 + Power(x,2)))*Sin(Pi*x))*Sin(4*Pi*y) +(45*Sin(Pi*t)*(1 +18*Power(Pi,2) -18*Power(Pi,2)*Power(x,2) +(-1 +26*Power(Pi,2)*(-1 + Power(x,2)))*Cos(4*Pi*x) +8*Pi*x*Sin(4*Pi*x))*Sin(6*Pi*y))/2.))/R+ (8*Sin((Pi*t)/2.)*Sin(Pi*x)*(Pi*(-1 + Power(x,2))*Cos(Pi*x) +x*Sin(Pi*x))*Sin(4*Pi*y) +45*Sin(Pi*t)*Sin(2*Pi*x)*(2*Pi*(-1 + Power(x,2))*Cos(2*Pi*x) +x*Sin(2*Pi*x))*Sin(6*Pi*y))*(8*Pi*(-1 + Power(x,2))*Sin((Pi*t)/2.)*Power(Sin(Pi*x),2)*Sin(4*Pi*y) +5*(-40 +9*Pi*(-1 + Power(x,2))*Sin(Pi*t)*Power(Sin(2*Pi*x),2)*Sin(6*Pi*y)))))/20000.,(400*Pi*(2*Cos((Pi*t)/2.)*Sin(Pi*x)*(Pi*(-1 + Power(x,2))*Cos(Pi*x) +x*Sin(Pi*x))*Power(Sin(2*Pi*y),2) +15*Cos(Pi*t)*Sin(2*Pi*x)*(2*Pi*(-1 + Power(x,2))*Cos(2*Pi*x) +x*Sin(2*Pi*x))*Power(Sin(3*Pi*y),2)) -(800*Pi*(16*Pi*Power(Cos(2*Pi*y),2)*Sin((Pi*t)/2.)*Sin(Pi*x)*(Pi*(-1 + Power(x,2))*Cos(Pi*x) +x*Sin(Pi*x)) +135*Pi*Power(Cos(3*Pi*y),2)*Sin(Pi*t)*Sin(2*Pi*x)*(2*Pi*(-1 + Power(x,2))*Cos(2*Pi*x) +x*Sin(2*Pi*x)) +12*Pi*x*Power(Cos(Pi*x),2)*Sin((Pi*t)/2.)*Power(Sin(2*Pi*y),2)- 28*Pi*x*Sin((Pi*t)/2.)*Power(Sin(Pi*x),2)*Power(Sin(2*Pi*y),2)+ 6*Sin((Pi*t)/2.)*Sin(2*Pi*x)*Power(Sin(2*Pi*y),2)+ 12*Power(Pi,2)*Sin((Pi*t)/2.)*Sin(2*Pi*x)*Power(Sin(2*Pi*y),2)- 12*Power(Pi,2)*Power(x,2)*Sin((Pi*t)/2.)*Sin(2*Pi*x)*Power(Sin(2*Pi*y),2)+ 180*Pi*x*Power(Cos(2*Pi*x),2)*Sin(Pi*t)*Power(Sin(3*Pi*y),2)- 315*Pi*x*Sin(Pi*t)*Power(Sin(2*Pi*x),2)*Power(Sin(3*Pi*y),2)+ 45*Sin(Pi*t)*Sin(4*Pi*x)*Power(Sin(3*Pi*y),2)+ 255*Power(Pi,2)*Sin(Pi*t)*Sin(4*Pi*x)*Power(Sin(3*Pi*y),2)- 255*Power(Pi,2)*Power(x,2)*Sin(Pi*t)*Sin(4*Pi*x)*Power(Sin(3*Pi*y),2)))/R + 4*Pi*(4*Sin((Pi*t)/2.)*Sin(Pi*x)*(Pi*(-1 + Power(x,2))*Cos(Pi*x) +x*Sin(Pi*x))*Power(Sin(2*Pi*y),2) +15*Sin(Pi*t)*Sin(2*Pi*x)*(2*Pi*(-1 + Power(x,2))*Cos(2*Pi*x) +x*Sin(2*Pi*x))*Power(Sin(3*Pi*y),2))*(8*Sin((Pi*t)/2.)*Sin(Pi*x)*(Pi*(-1 + Power(x,2))*Cos(Pi*x) +x*Sin(Pi*x))*Sin(4*Pi*y) +45*Sin(Pi*t)*Sin(2*Pi*x)*(2*Pi*(-1 + Power(x,2))*Cos(2*Pi*x) +x*Sin(2*Pi*x))*Sin(6*Pi*y)) -(8*Power(Pi,2)*(-1 + Power(x,2))*Power(Cos(Pi*x),2)*Sin((Pi*t)/2.)*Power(Sin(2*Pi*y),2) -8*Sin((Pi*t)/2.)*Sin(Pi*x)*(-4*Pi*x*Cos(Pi*x) +(-1 +Power(Pi,2)*(-1 + Power(x,2)))*Sin(Pi*x))*Power(Sin(2*Pi*y),2) +15*Sin(Pi*t)*(1 +(-1 +8*Power(Pi,2)*(-1 + Power(x,2)))*Cos(4*Pi*x) +8*Pi*x*Sin(4*Pi*x))*Power(Sin(3*Pi*y),2))*(8*Pi*(-1 + Power(x,2))*Sin((Pi*t)/2.)*Power(Sin(Pi*x),2)*Sin(4*Pi*y) +5*(-40 +9*Pi*(-1 + Power(x,2))*Sin(Pi*t)*Power(Sin(2*Pi*x),2)*Sin(6*Pi*y))))/40000.);

      // Define our solution
      Expr uExact = List(1 - (Pi*(-1 + Power(x,2))*(8*Sin((Pi*t)/2.)*Power(Sin(Pi*x),2)*Sin(4*Pi*y) +45*Sin(Pi*t)*Power(Sin(2*Pi*x),2)*Sin(6*Pi*y)))/200.,(4*Sin((Pi*t)/2.)*Sin(Pi*x)*(Pi*(-1 + Power(x,2))*Cos(Pi*x) + x*Sin(Pi*x))*Power(Sin(2*Pi*y),2) +15*Sin(Pi*t)*Sin(2*Pi*x)*(2*Pi*(-1 + Power(x,2))*Cos(2*Pi*x) + x*Sin(2*Pi*x))*Power(Sin(3*Pi*y),2))/100.);

      Expr pExact = 0.0;


      

      // Define our mesh
      // Might be a good idea to read this in from the ForwardProblem 
      MeshType spatialMeshType = new Sundance::BasicSimplicialMeshType();
      double xmin = 0.0;
      double xmax = 1.0;
      double ymin = 0.0;
      double ymax = 1.0;
      MeshSource spatialMesher = new Sundance::PartitionedRectangleMesher(xmin,xmax,nx,1,ymin,ymax,nx,1,spatialMeshType,0);
      Mesh spatialMesh = spatialMesher.getMesh();



      // Specify the location of the velocity snapshots
      string outDir = "/home/sirush/PhDResearch/MMS_Transient_Channel/ForwardProblem/Results/Re";
      std::ostringstream ReynoldsString;
      ReynoldsString << std::setprecision(precision) << Re;
      outDir += ReynoldsString.str() + "/nx" +Teuchos::toString(nx) + "nt" + Teuchos::toString(nSteps) + "/";
      //      string tag = "st-v";
      string filePrefix = outDir + "st-v";

      // Specify the location of the reduced-order basis functions
      string POD_DataDir = "/home/sirush/PhDResearch/MMS_Transient_Channel/POD/SingleParameterSpace/Results/Re" + ReynoldsString.str() + "/nx" + Teuchos::toString(nx) + "nt" + Teuchos::toString(nSteps) + "/tol";
      std::ostringstream tolFileValue;
      tolFileValue << std::setprecision(precision) << tol;
      POD_DataDir = POD_DataDir + tolFileValue.str() + "/";

      // Create a BasisFamily to express our solution; each spatial dim is Lagrange(2)
      Array<Sundance::BasisFamily> velBasis(spatialMesh.spatialDim());
      for(int i = 0; i < velBasis.length(); i++)
	velBasis[i] = new Sundance::Lagrange(2);
      
      // Define our vector type
      Playa::VectorType<double> epetraVecType = new Playa::EpetraVectorType();
      Sundance::DiscreteSpace velocityDS(spatialMesh, velBasis, epetraVecType);

      string NLParamFile = "playa-newton-armijo.xml";

      // Calculate the kinematic viscosity nu
      double nu = (U0*L0)/Re;
      
      // Create a velocityROM object
      velocityROM ROM(POD_DataDir,NLParamFile,velocityDS, u0, q,
		      t, nSteps, deltat, nu, tol, verbosity);
      ROM.initialize(filePrefix);
      ROM.generate_alpha();
      // Get the discretized in time Array of velocity functions
      Array<Expr> uRO(ROM.get_uRO() );

      // Write alphaROM(t) to file
      string ROM_Dir="Results/SingleParameterSpace/tol" + tolFileValue.str()
	+ "/Re" + ReynoldsString.str() + "/nx"+Teuchos::toString(nx)+"nt"
	+Teuchos::toString(nSteps)+"/";
      string alpha_filename = ROM_Dir + "alphaROM.txt";
      int dirCreation = system( ("mkdir -p " + ROM_Dir).c_str() );
      TEUCHOS_TEST_FOR_EXCEPTION( dirCreation == -1, runtime_error,
				  "Failed to create " + ROM_Dir);
      ROM.write_alpha(alpha_filename);


	  

      // Find the error between alphaEx and alphaROM
      Vector<double> alphaError = ROM.alphaErrorCheck(uExact);
      cout << "Run for nx = " << nx << ", nSteps = " << nSteps << endl;
      cout << "||2norm of the error in alpha at all timesteps||_2  :\t "  << alphaError.norm2() << endl;
      cout << "||2norm of the error in alpha at all timesteps||_inf:\t " << alphaError.normInf() << endl;
      

      Vector<double> uError = ROM.velocityErrorCheck(uExact);

      SUNDANCE_ROOT_MSG1(verbosity, "Number of velocity modes kept: " + Teuchos::toString(ROM.get_phi().size()));
      Out::root() << "||2norm of the error in u at all timesteps||_2  :\t " << uError.norm2() << endl;
      Out::root() << "||2norm of the error in u at all timesteps||_inf:\t " << uError.normInf() << endl;


      timer.stop();
      Out::root() << "runtime=" << timer.totalElapsedTime() << endl << endl;
      
      // Visualize the results
      SUNDANCE_ROOT_MSG1(verbosity, "Writing results to file");
      string vtkDir = ROM_Dir;
      string vtkfilename = "Re" + ReynoldsString.str() + "nx"+Teuchos::toString(nx)+"nt"+Teuchos::toString(nSteps);
      //      vtkDir = vtkDir + vtkfilename + "/";


      DiscreteSpace scalarDS(spatialMesh, new Lagrange(1), new EpetraVectorType());
      L2Projector projector(velocityDS, uExact);
      //Write uExact for all the time steps
      for(int time=0; time< nSteps+1; time++)
	{
	  FieldWriter writer = new VTKWriter(vtkDir+vtkfilename+"step"+Teuchos::toString(time));
	  writer.addMesh(spatialMesh);
	  t.setParameterValue(tInit+time*deltat);
	  L2Projector projectorRO(velocityDS, uRO[time]);
	  L2Projector uErrorProjector(velocityDS, uExact - uRO[time]);
	  Expr absErr = sqrt( (uExact - uRO[time])*(uExact - uRO[time]));
	  Expr absU = sqrt(uExact * uExact);
	  Expr absURO = sqrt(uRO[time] * uRO[time]);
	  L2Projector uMagProj(scalarDS, absU);
	  L2Projector uROMagProj(scalarDS, absURO);
	  L2Projector absErrorProj(scalarDS, absErr);
	  L2Projector relErrorProj(scalarDS, absErr / (absU + 1.0));
	  writer.addField("uMag", new ExprFieldWrapper(uMagProj.project()[0]) );
	  writer.addField("uROMag", new ExprFieldWrapper(uROMagProj.project()[0]) );
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

