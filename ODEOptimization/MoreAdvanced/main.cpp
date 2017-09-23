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
      Time timer("total");
      timer.start();
      
      
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

      // int Ru = 2;
      // Sundance::setOption("Ru", Ru, "Number of POD basis functions");
      double tol = 0.999;
      Sundance::setOption("tol",tol,"The tolerance used in the RIC to select the number of reduced-order basis functions to keep");

      double tFinal = 1.0;
      Sundance::setOption("tFinal",tFinal,"Final time value");

      int meshVerb = 0;
      Sundance::setOption("meshVerb",meshVerb,"Mesh verbosity level");

      double eta = 0.05;
      Sundance::setOption("eta",eta,"constant parameter on the regularization term");

      double eta_design = 0.001;
      Sundance::setOption("eta_design",eta_design,"Value for the constant term in the design equation");
      
      double eta_reg = 10.0;
      Sundance::setOption("eta_reg",eta_reg,"Value for the constant term in the regularization term");

      Sundance::init(&argc, &argv);

      // State the location of the directory holds alphaROM
      string ROM_base_dir = "/home/sirush/PhDResearch/MMS_Transient_Channel/Results/ROM/uRO/nx" + Teuchos::toString(nx) + "nt" + Teuchos::toString(nSteps) + "/";

      // State the location of the POD files (phi, b, A, and T)
      string POD_DataDir = "/home/sirush/PhDResearch/MMS_Transient_Channel/Results/POD/nx" + Teuchos::toString(nx) + "nt" + Teuchos::toString(nSteps) + "/tol" + std::to_string(tol) + "/";

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
      string snapshotDataDir = "/home/sirush/PhDResearch/MMS_Transient_Channel/Results/ForwardProblem/R=1forward_problem_TransientChannel_nx" + Teuchos::toString(nx) + "-nt-" + Teuchos::toString(nSteps) + "/";
      string filenamePrefix = "st-v";
      string snapshotFilenamePrefix = snapshotDataDir + filenamePrefix;

      // Create the sensorData object that will generate vstar
      sensorData sensorObj(positionArray, spatialMesh, timeMesh, snapshotFilenamePrefix, nSteps, eastVec, quadOrder);
     

      // Define the KKT object for this test problem (Transient Channel)
      KKT_Transient_Channel KKT_System(POD_DataDir, spatialMesh, timeMesh, sensorObj, tFinal, nSteps, quadOrder, verbosity);
      KKT_System.initialize();

      // Solve the KKT system
      string solverFile = "playa-newton-amesos.xml";
      Expr soln = KKT_System.solve(solverFile,eta_design,eta_reg);

      Expr alphaOPT = soln[0];
      for(int r=1; r < 2; r++)
      	alphaOPT.append(soln[r]);

      double alphaError = KKT_System.errorCheck(ROM_base_dir, alphaOPT);
      cout << "Run for nx = " << nx << ", nSteps = " << nSteps << ", and tol = " << tol << endl;
      cout << "||alphaROM - alphaOPT||_2 = " << alphaError << endl;



      timer.stop();
      Out::root() << "runtime = " << timer.totalElapsedTime() << endl << endl;

      /* The solve for alpha values will be the first Ru components of U0
       We need to now build the approximation to the velocity uOpt
      Expr uOpt = uB;
      for(int r = 0; r < Ru; r++)
	{
	  uOpt = uOpt + U0[r]*phi[r];
	}
      */
      
      
      // FieldWriter writer = new DSVWriter("2DLinearOpt-.dat");
      // writer.addMesh(mesh);
      // writer.addField("x[0]", new ExprFieldWrapper(U0[0]));
      // writer.addField("x[1]", new ExprFieldWrapper(U0[1]));
      // writer.addField("lambda[0]", new ExprFieldWrapper(U0[2]));
      // writer.addField("lambda[1]", new ExprFieldWrapper(U0[3]));
      // writer.addField("alpha[0]", new ExprFieldWrapper(U0[4]));
      // writer.addField("alpha[1]", new ExprFieldWrapper(U0[5]));
      // writer.write();
      
      // Array<double> alphaNum = Teuchos::tuple(
      // 					      L2Norm(mesh, left, U0[0], quad),
      // 					      L2Norm(mesh, left, U0[1], quad)
      // 					      );
      // for (int j=0; j<2; j++)
      // 	{
      // 	  Tabs tab1;
      // 	  Out::os() << tab1 << "Alpha[" << j << "]: exact="
      // 		    << alphaExact[j]
      // 		    << ", numerical " << alphaNum[j]
      // 		    << ", error=" << fabs(alphaExact[j] - alphaNum[j])
      // 		    << endl;
      // 	}


     
    }
  catch(std::exception& ex)
    {
      cerr << "main() caught exception: " << ex.what() << endl;
    }
  Sundance::finalize();
}

    

    
    





