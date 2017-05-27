//Local files
#include "PlayaSVD.hpp"
// Kathryn, this is included so that Pi has a definition
#include "MathematicaConverter.hpp"

//Viento files
#include "VientoSnapshotIO.hpp"

using std::cout;
using std::endl;

int main(int argc, char *argv[]) 
{
  try
    {
      Time timer("total");
      timer.start();
      
      int verbosity = 1;	
      Sundance::setOption("verbosity", verbosity, "verbosity level");

      int nx = 16;
      Sundance::setOption("nx", nx, "Number of elements along each axis");

      int nSteps = 16;
      Sundance::setOption("nSteps", nSteps, "Number of time steps");

      Sundance::init(&argc, &argv);

      // Define our mesh
      MeshType meshType = new Sundance::BasicSimplicialMeshType();
      double xmin = 0.0;
      double xmax = 2.0*Pi;
      double ymin = 0.0;
      double ymax = 2.0*Pi;
      MeshSource mesher = new Sundance::PartitionedRectangleMesher(xmin,xmax,nx,1,ymin,ymax,nx,1,meshType,0);
      Mesh mesh = mesher.getMesh();



      // Read the snapshots into a matrix
      // Kathryn, you will probably need to change the file directory for the results of the homogeneous Dirichlet problem
      string fileDir = "Results/ForwardProblem/HG_Dirichlet_nx" + Teuchos::toString(nx) + "_nt" + Teuchos::toString(nSteps);
      string tag = "st-v";
      string filename = fileDir + "/" + tag;

      // Create a BasisFamily to express our solution
      Array<Sundance::BasisFamily> ubasis = List(new Sundance::Lagrange(2), new Sundance::Lagrange(2)); // 2nd order Piece-Wise Quad Lagrange in 2D
      // Define our vector type
      Playa::VectorType<double> epetraVecType = new Playa::EpetraVectorType();
      Sundance::DiscreteSpace uDS(mesh, ubasis, epetraVecType);

      // Perform the POD of the matrix W
      
      LinearOperator<double> W = snapshotToMatrix(filename, nSteps, mesh);
      Playa::LinearOperator<double> U;
      Playa::LinearOperator<double> Phi;
      Playa::Vector<double> sigma;


      // W and mesh need to be defined
      SUNDANCE_ROOT_MSG1(verbosity, "Entering POD for velocity");
      POD(W,sigma,U,Phi,uDS,verbosity);
      SUNDANCE_ROOT_MSG2(verbosity, "POD finished for velocity");

      
      timer.stop();
      cout << "runtime=" << timer.totalElapsedTime() << endl;


    }
  catch(std::exception& e)
    {
      Out::root() << "Caught exception: " << e.what() << std::endl;
      return -1;
    }
}
