#include "Sundance.hpp"
#include "VientoDefaultOutputManager.hpp"
#include "VientoSnapshotIO.hpp"
#include "PDEOptPointData.hpp"

// Standard Library Functions
using std::cout;
using std::endl;
using std::vector;

/**
 * The purpose of this program is to test how to create a PointData object and use it 
 * to depict our Lidar sensors. Note: the test functions are symbolic, and thus they must first be
 * projected in order to turn them into DiscreteFunciton objects.
 */

/**
 * getPoint(Point P) returns a CellFilter that will return the cell containing
 * the P(x,y) point
 */
CellFilter getPoint(Point P)
{
  CellFilter vertices = new DimensionalCellFilter(0);
  CellFilter xPos = vertices.coordSubset(0, P[0]);
  CellFilter yPos = vertices.coordSubset(1, P[1]);
  CellFilter vertex = xPos.intersection(yPos);

  return vertex;
}

int main(int argc, char *argv[])
{
  try
    {
      int nx = 25;
      Sundance::setOption("nx",nx,"Number of spatial elements along the x-axis");

      int nSteps = 25;
      Sundance::setOption("nSteps", nSteps, "Number of timesteps taken to get to T");

      double tol;
      Sundance::setOption("tol", tol, "Value for the tolerance in the SloppyPointComparitor");
      
      
      Sundance::init(&argc, &argv);

      // Create the coordinate variables
      Expr x = new CoordExpr(0,"x");
      Expr y = new CoordExpr(1,"y");

      // Create the expression we wish to write to file
      Expr test1 = x+y;

      // Create the discrete space
      VectorType<double> vecType = new EpetraVectorType();
      MeshType meshType = new BasicSimplicialMeshType();
      double xmin = 0.0;
      double xmax = 1.0;
      double ymin = 0.0;
      double ymax = 1.0;
      MeshSource mesher = new Sundance::PartitionedRectangleMesher(xmin,xmax,nx,1,ymin,ymax,nx,1,meshType,0);
      Mesh mesh = mesher.getMesh();

      // Create an array that holds the * position values of the sensors
      Array<double> positionValues = tuple(0.2, 0.4, 0.6, 0.8);
      Array<Point> positionArray;
      positionArray.resize(pow(positionValues.length(),2));
      int numberOfReadings = positionArray.length();
      for(int i=0; i<positionValues.length(); i++)
	for(int j=0; j<positionValues.length(); j++)
	  {
	    positionArray[j+positionValues.length()*i].resize(2);
	    positionArray[j+positionValues.length()*i][0] = positionValues[i];
	    positionArray[j+positionValues.length()*i][1] = positionValues[j];
	  }
      
      
      cout << "Size of positionArray: " << positionArray.length() << endl;
      cout << "Values in positionArray: " << positionArray << endl << endl;

      // Create the east unit vector
      double eastAngle = 0.0;
      Expr eastVec = List(cos(eastAngle), sin(eastAngle));

      // Comeback and read a snapshot into an Expr, then dot it with eastVec
      string snapshotDataDir = "/home/sirush/PhDResearch/MMS_Transient_Channel/Results/ForwardProblem/R=1forward_problem_TransientChannel_nx" + Teuchos::toString(nx) + "-nt-" + Teuchos::toString(nSteps) + "/";
      string filenamePrefix = "st-v";
      string snapshotFilename = snapshotDataDir + filenamePrefix;
      int timeStep = 0;
      Expr velocityDF = readSnap(snapshotFilename, timeStep, mesh);

      // Calculate the integral in the east direction to establish the Array of sensor values
      QuadratureFamily quad = new GaussianQuadrature(6);
      CellFilter tempFilter;
      Array<double> valueArray;
      valueArray.resize(numberOfReadings);

      /* This works for taking the values from the forward simulation and storing them as the measurement resutls 
      for(int i = 0; i < numberOfReadings; i++)
	{
	  tempFilter = getPoint(positionArray[i]);
	  FunctionalEvaluator testIntegral(mesh, Integral(tempFilter, eastVec*velocityDF, quad));
	  valueArray[i] = testIntegral.evaluate();
	  cout << "The value of the integral eastVec*velocityDF over the point " << positionArray[i] << " is " << valueArray[i] << endl;
	}
	*/

      /* This works for setting the i^{th} mesasurement value to i */
      for(int i = 0; i < numberOfReadings; i++)
	{
	  valueArray[i] = i;
	  cout << "The value of the integral eastVec*velocityDF over the point " << positionArray[i] << " is " << valueArray[i] << endl;
	}
      

      // Create the PointData object
      PointData ptData(positionArray, valueArray, tol);
      CellFilter sensors = ptData.sensorLocations();
      Expr vstar = ptData.sensorValues();

      cout << "What vstar displays: " << endl << vstar << endl;
      for(int i = 0; i < numberOfReadings; i++)
	{
	  FunctionalEvaluator integral(mesh, Integral(getPoint(positionArray[i]), vstar, quad));
	  cout << "Attempting to access the measurement value for " << positionArray[i] << ": " << integral.evaluate() << endl;
	}

      // Integrating over all the measurement locations sums together their values
      FunctionalEvaluator all(mesh, Integral(sensors, vstar, quad));
	cout << "What happens if I integrate over all the measurement locations: " << all.evaluate() << endl;
      


    }
  catch(std::exception& e)
    {
      cerr << "main() caught exception: " << e.what() << endl;
    }
  Sundance::finalize();  
}
