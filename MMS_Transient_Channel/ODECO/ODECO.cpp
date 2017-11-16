#include "Sundance.hpp"
//#include "PlayaAmesosSolver.hpp"
#include "PlayaDenseSerialMatrix.hpp"
//#include "PlayaSerialVectorType.hpp"
#include "VientoSnapshotIO.hpp"
#include "denseSerialMatrixIO.hpp"


#include "integralOperator.hpp"
//#include "ODECO.hpp"
#include "KKT_Transient_Channel.hpp"
//#include "sensorData.hpp"

// Local Files
#include "MathematicaConverter.hpp"


// Standard Library Functions
using std::cout;
using std::endl;
using std::vector;

// Needed to write only significant digits of tol as part of the file name
#include <sstream>
#include <iomanip>

/**
 * getPoint(Point P) returns a CellFilter that will return the cell containing
 * the P(x,y) point
 */
// CellFilter getPoint(Point P)
// {
//   CellFilter vertices = new DimensionalCellFilter(0);
//   CellFilter xPos = vertices.coordSubset(0, P[0]);
//   CellFilter yPos = vertices.coordSubset(1, P[1]);
//   CellFilter vertex = xPos.intersection(yPos);

//   return vertex;
// }


// Resume by cleaning up this file; and get results to be outputted.

    
int main(int argc, char *argv[]) 
{
  try
    {
      Time timerTotal("total");
      timerTotal.start();

      Time timerKKT("KKT");
      timerKKT.start();
      
      
      int verbosity = 0;      
      Sundance::setOption("verbosity", verbosity, "verbosity level");
      
      int nx = 25;
      Sundance::setOption("nx", nx, "Number of elements along each axis");

      int nSteps = 25;
      Sundance::setOption("nSteps", nSteps, "Number of time steps");

      //double gamma = 0.1;
      //Sundance::setOption("gamma", gamma, "perturbation to target");

      int quadOrder = 2;
      Sundance::setOption("quadOrder", quadOrder, "Order for the Gaussian Quadrature rule");

      double tol = 0.999;
      Sundance::setOption("tol",tol,"The tolerance used in the RIC to select the number of reduced-order basis functions to keep");

      int precision = 3;
      Sundance::setOption("precision",precision,"Number of significant digits to keep in a filename using tol");

      double Re = 1.0;
      Sundance::setOption("Re", Re, "Reynolds number");     

      double tFinal = 1.0;
      Sundance::setOption("tFinal",tFinal,"Final time value");

      int meshVerb = 0;
      Sundance::setOption("meshVerb",meshVerb,"Mesh verbosity level");

      double eta_design = 10.0;
      Sundance::setOption("eta_design",eta_design,"Value for the constant term in the design equation");
      
      double eta_reg = 0.001;
      Sundance::setOption("eta_reg",eta_reg,"Value for the constant term in the regularization term");

      string parameterSpace = "single";
      Sundance::setOption("parameterSpace",parameterSpace,"Defines the space from which the POD velocity modes are derived; i.e. a single selection of parameters or multiple simulations");

      Sundance::init(&argc, &argv);

      // State the location of the POD files (phi, b, A, and T)
      std::ostringstream ReynoldsString;
      ReynoldsString << std::setprecision(precision) << Re;
      string POD_DataDir = "../POD/";
      if(parameterSpace=="single")
	{
	  cout << "Took single branch" << endl;
	  POD_DataDir += "SingleParameterSpace/";
	}
      else if(parameterSpace=="multiple")
	{
	  cout << "Took multiple branch" << endl;
	}
      POD_DataDir += "Results/Re" + ReynoldsString.str()
	+ "/nx" + Teuchos::toString(nx) + "nt" + Teuchos::toString(nSteps)
	+ "/tol";
      std::ostringstream tolFileValue;
      tolFileValue << std::setprecision(precision) << tol;
      POD_DataDir += tolFileValue.str() + "/";

      // Create the spatial mesh
      MeshType spatialMeshType = new Sundance::BasicSimplicialMeshType();
      double xmin = 0.0;
      double xmax = 1.0;
      double ymin = 0.0;
      double ymax = 1.0;
      MeshSource spatialMesher = new Sundance::PartitionedRectangleMesher(xmin,xmax,nx,1,ymin,ymax,nx,1,spatialMeshType,meshVerb);
      Mesh spatialMesh = spatialMesher.getMesh();

      // Define the time mesh
      MeshType timeMeshType = new BasicSimplicialMeshType();
      MeshSource timeMesher = new PartitionedLineMesher(0.0, tFinal, nSteps, timeMeshType);
      Mesh timeMesh = timeMesher.getMesh();



      // Create the measurement value functions
      // Create an array that holds the * position values of the sensors
      Array<double> positionValues = tuple(0.2, 0.4, 0.6, 0.8);
      Array<Point> positionArray;
      // spatialDim() returns n for nD
      positionArray.resize(pow(positionValues.length(),spatialMesh.spatialDim() ));
      // Assigns values for (x1,y1) through (x1,yN), then goes to x2 and so forth
      for(int i=0; i<positionValues.length(); i++)
	for(int j=0; j<positionValues.length(); j++)
	  {
	    positionArray[j+positionValues.length()*i].resize(spatialMesh.spatialDim() );
	    positionArray[j+positionValues.length()*i][0] = positionValues[i];
	    positionArray[j+positionValues.length()*i][1] = positionValues[j];	  
	  }

      // Create the sample direciton vector
      double eastAngle = 0.0;
      Expr eastVec = List(cos(eastAngle), sin(eastAngle));

      // Specify where to find the sensor data for one choice of parameter values
      string snapshotDataDir = "../ForwardProblem/Results/Re"
	+ ReynoldsString.str() + Teuchos::toString(nx) + "nt" + Teuchos::toString(nSteps) + "/";
      string filenamePrefix = "st-v";
      string snapshotFilenamePrefix = snapshotDataDir + filenamePrefix;

      // Create the sensorData object that will generate vstar
      sensorData sensorObj(positionArray, spatialMesh, timeMesh, snapshotFilenamePrefix, nSteps, eastVec, quadOrder);
     

      // Define the KKT object for this test problem (Transient Channel)
      KKT_Transient_Channel KKT_System(POD_DataDir, spatialMesh, timeMesh,
				       sensorObj, tFinal, nSteps, quadOrder, verbosity);
      KKT_System.initialize();

      // Solve the KKT system
      string solverFile = "playa-newton-amesos.xml";
      Expr soln = KKT_System.solve(solverFile,eta_design,eta_reg);

      Expr alphaOPT  = KKT_System.get_alpha();
      Array<Expr> uOPT = KKT_System.get_uOPT();

      timerKKT.stop();


      // Check alphaOPT against alphaExact
      Time timerErrorCheck("KKT");
      timerErrorCheck.start();
      Expr x = new CoordExpr(0,"x");
      Expr y = new CoordExpr(1,"y");
      Expr t = new Sundance::Parameter(0.0);
      Expr uExact = List(1 - (Pi*(-1 + Power(x,2))*(8*Sin((Pi*t)/2.)*Power(Sin(Pi*x),2)*Sin(4*Pi*y) +45*Sin(Pi*t)*Power(Sin(2*Pi*x),2)*Sin(6*Pi*y)))/200.,(4*Sin((Pi*t)/2.)*Sin(Pi*x)*(Pi*(-1 + Power(x,2))*Cos(Pi*x) + x*Sin(Pi*x))*Power(Sin(2*Pi*y),2) +15*Sin(Pi*t)*Sin(2*Pi*x)*(2*Pi*(-1 + Power(x,2))*Cos(2*Pi*x) + x*Sin(2*Pi*x))*Power(Sin(3*Pi*y),2))/100.);

      // Returns [absolute alpha error, relative alpha error,
      //          aggregate velocity abs error, aggregate velocity rel error]
      Array<double> errors = KKT_System.errorCheck(uExact,t);
      cout << "Re=" << Re
	   << " nx=" << nx
	   << " nSteps=" << nSteps
	   << " tol=" << tol << endl;
      cout << "||alphaExact - alphaOPT||_2 = " << errors[0] << endl;
      cout << "||alphaExact - alphaOPT||_2 / ||alphaExact||_2 = " << errors[1] << endl;
      cout << "||L2norm(uEx(t_m) - uOPT_(t_m)) at all timesteps||_2 :\t "  << errors[2] << endl;
      cout << "||Relative error in velocity at all timesteps||_2 :\t " << errors[3] << endl;


      timerErrorCheck.stop();

      cout << "Writing results to file" << endl;
      Time timerVTK("vtk");
      timerVTK.start();
      /*
      string vtkDir = "Results/";
      if (parameterSpace=="single")
      	vtkDir+="SingleParameterSpace/";
      else if (parameterSpace=="multiple")
      	vtkDir+="MultipleParameterSpace/";

      vtkDir += "tol" + tolFileValue.str() + "/Re" + ReynoldsString.str() + "/"
      	+ "nx" + Teuchos::toString(nx) + "nt" + Teuchos::toString(nSteps) + "/";

      // Create the directory for storing the vtk files
      int dirCreation = system( ("mkdir -p " + vtkDir).c_str() );
      TEUCHOS_TEST_FOR_EXCEPTION( dirCreation == -1,
      				  runtime_error,
      				  "Failed to create " + vtkDir);
      
	    
      DiscreteSpace scalarDS(spatialMesh, new Lagrange(1), new EpetraVectorType());
      // Create a BasisFamily to express our solution; each spatial dim is Lagrange(2)
      Array<Sundance::BasisFamily> velBasis(spatialMesh.spatialDim());
      for(int i = 0; i < velBasis.length(); i++)
      	velBasis[i] = new Sundance::Lagrange(2);
      DiscreteSpace velocityDS(spatialMesh, velBasis, new EpetraVectorType());
      L2Projector projector(velocityDS, uExact);
      //Write uExact for all the time steps
      int tInit = 0;
      for(int time=0; time< nSteps+1; time++)
      	{
      	  FieldWriter writer = new VTKWriter(vtkDir+"step"+Teuchos::toString(time));
      	  writer.addMesh(spatialMesh);
      	  t.setParameterValue(tInit+time*(tFinal-tInit)/nSteps);
      	  L2Projector projectorOPT(velocityDS, uOPT[time]);
      	  L2Projector uErrorProjector(velocityDS, uExact - uOPT[time]);
      	  Expr absErr = sqrt( (uExact - uOPT[time])*(uExact - uOPT[time]));
      	  Expr absU = sqrt(uExact * uExact);
      	  Expr absUOPT = sqrt(uOPT[time] * uOPT[time]);
      	  L2Projector uMagProj(scalarDS, absU);
      	  L2Projector uOPTMagProj(scalarDS, absUOPT);
      	  L2Projector absErrorProj(scalarDS, absErr);
      	  L2Projector relErrorProj(scalarDS, absErr / absU);
      	  writer.addField("uMag", new ExprFieldWrapper(uMagProj.project()[0]) );
      	  writer.addField("uOPTMag", new ExprFieldWrapper(uOPTMagProj.project()[0]) );
      	  writer.addField("errAbs", new ExprFieldWrapper(absErrorProj.project()[0]) );
      	  writer.addField("errRel", new ExprFieldWrapper(relErrorProj.project()[0]) );
      	  writer.addField("uExact[0]", new ExprFieldWrapper(projector.project()[0]) );
      	  writer.addField("uExact[1]", new ExprFieldWrapper(projector.project()[1]) );
      	  writer.addField("uOPT[0]", new ExprFieldWrapper(projectorOPT.project()[0]) );
      	  writer.addField("uOPT[1]", new ExprFieldWrapper(projectorOPT.project()[1]) );
      	  writer.addField("uError[0]", new ExprFieldWrapper(uErrorProjector.project()[0]) );
      	  writer.addField("uError[1]", new ExprFieldWrapper(uErrorProjector.project()[1]) );
      	  writer.write();	  
      	}
      */
      timerVTK.stop();
      timerTotal.stop();
      Out::root() << "KKT Runtime = " << timerKKT.totalElapsedTime() << endl;
      Out::root() << "ErrorCheck Runtime = " << timerErrorCheck.totalElapsedTime() << endl;
      Out::root() << "VTK Runtime = " << timerVTK.totalElapsedTime() << endl;
      Out::root() << "Total Runtime = " << timerTotal.totalElapsedTime() << endl;
     
    }
  catch(std::exception& ex)
    {
      cerr << "main() caught exception: " << ex.what() << endl;
    }
  Sundance::finalize();
}

    

    
    





