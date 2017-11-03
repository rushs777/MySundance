#include "Sundance.hpp"
#include "POD_SVD.hpp" // For the POD
#include "VientoSnapshotIO.hpp" // for snapshotToMatrix
#include "denseSerialMatrixIO.hpp" // To get uB

// Needed to write only significant digits of tol as part of the file name
#include <sstream>
#include <iomanip>

using std::cout;
using std::endl;

int main(int argc, char* argv[])
{
  Time timer("total");
  timer.start();
  
  int nx = 8;
  Sundance::setOption("nx",nx,"Number of spatial elements along each axis");

  int nSteps = 8;
  Sundance::setOption("nSteps",nSteps,"Number of time steps taken to go from iInit to tFinal");

  double tol = .999;
  Sundance::setOption("tol",tol,"The tolerance used in the relative information content (RIC) to determine the number of basis functions to keep");

  int precision = 3;
  Sundance::setOption("precision",precision,"Number of significant digits to keep in a filename using tol");

  double Re = 1.0;
  Sundance::setOption("Re", Re, "Reynolds number");

  int verbosity = 1;
  Sundance::setOption("verbosity",verbosity,"verbosity sets the level of displayed output");

  // Initialize Sundance
  Sundance::init(&argc, &argv);

  // Define the location of the snapshot matrix (matrices)
  string snapshotDataDir = "/home/sirush/PhDResearch/MMS_Transient_Channel/Results/ForwardProblem/Re";
  std::ostringstream ReynoldsString;
  ReynoldsString << std::setprecision(precision) << Re;
  snapshotDataDir = snapshotDataDir + ReynoldsString.str()
    + "/nx" + Teuchos::toString(nx) + "nt" + Teuchos::toString(nSteps) + "/";  

  
  string velocityPrefix = "st-v";
  string snapshotFilenamePrefix = snapshotDataDir + velocityPrefix;

  // Create the spatial mesh
  // Change this to reading in the mesh
  MeshType spatialMeshType = new Sundance::BasicSimplicialMeshType();
  double xmin = 0.0;
  double xmax = 1.0;
  double ymin = 0.0;
  double ymax = 1.0;
  MeshSource spatialMesher = new Sundance::PartitionedRectangleMesher(xmin,xmax,nx,1,
								      ymin,ymax,nx,1,
								      spatialMeshType,0);
  Mesh spatialMesh = spatialMesher.getMesh();  

  
  // Read in the matrix
  LinearOperator<double> W = snapshotToMatrix(snapshotFilenamePrefix, nSteps, spatialMesh);

  // Create the DiscreteSpace object. Lagrange(2) is used by the forward
  // problem simulation code for each spatial dimension
  Array<Sundance::BasisFamily> basisArray(spatialMesh.spatialDim());
  for(int i = 0; i < basisArray.length(); i++)
    basisArray[i] = new Sundance::Lagrange(2);

  VectorType<double> epetraVecType = new Playa::EpetraVectorType();
  DiscreteSpace velocityDS(spatialMesh, basisArray, epetraVecType);
  
  LinearOperator<double> Wprime = generateFluctuationMatrix(W);

  // Create the POD object
  POD_SVD pod(Wprime, velocityDS, verbosity);
  pod.calculateSVD();
  pod.calculateBasisFunctions();
  string POD_DataDir = "/home/sirush/PhDResearch/MMS_Transient_Channel/Results/POD/Re"
    + ReynoldsString.str() +"/nx" + Teuchos::toString(nx)
    + "nt" + Teuchos::toString(nSteps) + "/tol";
  std::ostringstream tolFileValue;
  tolFileValue << std::setprecision(precision) << tol;
  POD_DataDir = POD_DataDir + tolFileValue.str() + "/";
  Array<Expr> phi = pod.get_basis_functions(tol,POD_DataDir);

  timer.stop();
  cout << "For Re = " << Re << ", nx = " << nx << ", nSteps = " << nSteps << ", and tol = "
       << tol << endl;
  cout << "runtime = " << timer.totalElapsedTime() << endl << endl << endl;

  Sundance::finalize(); 

}
