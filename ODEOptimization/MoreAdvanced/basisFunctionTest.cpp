#include "Sundance.hpp"
#include "VientoDefaultOutputManager.hpp"
#include "VientoSnapshotIO.hpp"

// Standard Library Functions
using std::cout;
using std::endl;
using std::vector;

/**
 * This program tests using the VientoSnapshotIO functions to write DiscreteFunctions to file
 * and read them back. Note: the test functions are symbolic, and thus they must first be
 * projected in order to turn them into DiscreteFunciton objects.
 */

int main(int argc, char *argv[])
{
  try
    {
      int nx = 1;
      Sundance::setOption("nx",nx,"Number of spatial elements along the x-axis");
      
      
      Sundance::init(&argc, &argv);

      //RCP<DefaultOutputManager> output = rcp(new DefaultOutputManager(filename,

      // Create the coordinate variable
      Expr x = new CoordExpr(0,"x");

      // Create the expression we wish to write to file
      Expr test1 = x*x+3;

      // Create the discrete space
      VectorType<double> vecType = new EpetraVectorType();
      MeshType meshType = new BasicSimplicialMeshType();
      double xmin = 0.0;
      double xmax = 1.0;
      MeshSource mesher = new PartitionedLineMesher(xmin,xmax,nx,meshType,2);
      Mesh mesh = mesher.getMesh();

      BasisFamily bas = new Lagrange(2);
      DiscreteSpace ds(mesh,bas,vecType);
      L2Projector projector(ds,test1);
      Expr df_test1 = projector.project();
      Vector<double> test1Vec = getDiscreteFunctionVector(df_test1);

      cout << "test1: " << test1 << endl;
      cout << "test1 evaluated on the mesh: " << endl << df_test1 << endl;
      cout << "Vector of values: " << endl << test1Vec << endl;

      //DiscreteSpace ds = getDiscreteSpace(test1);
      string outDir = "Results/test/";
      string filename = outDir + "DF";
      int tag = 0;
      system( ("mkdir -p " + outDir).c_str() );
      writeSnap(filename,tag,df_test1);

      cout << "Attempting to read back test1" << endl;
      Expr test2 = readSnap(filename,tag,mesh);

      cout << "Value of test2: " << test2 << endl;

      cout << "getDiscreteFunctionVector(test2): " << endl << getDiscreteFunctionVector(test2) << endl;

      cout << "What happens if I say test1 - test2 ? " << endl << test1 - test2 << endl;

      CellFilter interior = new MaximalCellFilter();
      QuadratureFamily quad = new GaussianQuadrature(6);
      cout << "L2Norm of the difference: " << L2Norm(mesh, interior, test1-test2, quad) << endl;
      
    }
  catch(std::exception& e)
    {
      cerr << "main() caught exception: " << e.what() << endl;
    }
  Sundance::finalize();  
}
