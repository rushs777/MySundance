#include "Sundance.hpp"
#include "POD_SVD.hpp" // For the POD
#include "VientoSnapshotIO.hpp" // for snapshotToMatrix
#include "denseSerialMatrixIO.hpp" // To get uB, matrixAssembly

// Needed to write only significant digits of tol as part of the file name
#include <sstream>
#include <iomanip>

using std::cout;
using std::endl;

int main(int argc, char* argv[])
{
  Time timer("total");
  timer.start();
  
  int nx = 24;
  Sundance::setOption("nx",nx,"Number of spatial elements along each axis");

  int nSteps = 24;
  Sundance::setOption("nSteps",nSteps,"Number of time steps taken to go from iInit to tFinal");

  double tFinal = 1.0;
  Sundance::setOption("tFinal", tFinal, "final time");  

  double tol = .999;
  Sundance::setOption("tol",tol,"The tolerance used in the relative information content (RIC) to determine the number of basis functions to keep");

  int precision = 3;
  Sundance::setOption("precision",precision,"Number of significant digits to keep in a filename using tol");

  // double Re = 1.0;
  // Sundance::setOption("Re", Re, "Reynolds number");

  int verbosity = 1;
  Sundance::setOption("verbosity",verbosity,"verbosity sets the level of displayed output");

  // Initialize Sundance
  Sundance::init(&argc, &argv);

  // Specify the different Reynolds numbers to examine as an Array<double>
  // At this time, Sundance.hpp does not support passing Array<T> as command line arguments
  Array<int> ReValues = tuple(1,20,40,60,80,100);

  // Create the spatial mesh 
  // Change this to reading in the mesh (eventually)
  MeshType spatialMeshType = new Sundance::BasicSimplicialMeshType();
  double xmin = 0.0;
  double xmax = 1.0;
  double ymin = 0.0;
  double ymax = 1.0;
  MeshSource spatialMesher = new Sundance::PartitionedRectangleMesher(xmin,xmax,nx,1,
								      ymin,ymax,nx,1,
								      spatialMeshType,0);
  Mesh spatialMesh = spatialMesher.getMesh();

  // Create the DiscreteSpace object. Lagrange(2) is used by the forward
  // problem simulation code for each spatial dimension
  Array<Sundance::BasisFamily> basisArray(spatialMesh.spatialDim());
  for(int i = 0; i < basisArray.length(); i++)
    basisArray[i] = new Sundance::Lagrange(2);

  // Create the DiscreteSpace
  VectorType<double> epetraVecType = new Playa::EpetraVectorType();
  DiscreteSpace velocityDS(spatialMesh, basisArray, epetraVecType);


  /* 
   * Loop that will read in the snapshot matrices for a single choice of
   * simulation paramters (i.e. W_i)
   * Currently only accounts for different Reynolds number
   * Don't forget that nx has to be the same! (affects the number of rows)
   */

  // Variables needed from loop to loop
  // Define the location of the data from a single simulation (relative path)
  string snapshotDataDir;
  // Viento beings each velocity snapshot file with this prefix
  string velocityPrefix = "st-v";
  // The path to the velocity snapshot
  string snapshotFilenamePrefix;
  // Array object to hold the pointers to Wi, i=0:|Theta|
  Array<RCP<DenseSerialMatrix>> ptrLibrary;
  
  // Loop over Re
  for(int ReIndex = 0; ReIndex < ReValues.length(); ReIndex++)
    {
      snapshotDataDir = "../../ForwardProblem/Results/tFinal"
	+Teuchos::toString(int(tFinal)) +"sec/Re" + Teuchos::toString(ReValues[ReIndex])
	+ "/nx" + Teuchos::toString(nx) + "nt" + Teuchos::toString(nSteps) + "/";
      cout << "Reading in a result file from: " << endl << snapshotDataDir << endl;
      // Read in the matrix for a single simulation
      snapshotFilenamePrefix = snapshotDataDir + velocityPrefix;
      LinearOperator<double> Wi = snapshotToMatrix(snapshotFilenamePrefix, nSteps, spatialMesh);
      // Store a ptr to Wi
      ptrLibrary.push_back( DenseSerialMatrix::getConcretePtr(Wi) );
    }

  // From denseSerialMatrixIO.hpp use matrixAssembly to create a single matrix
  LinearOperator<double> W = matrixAssembly(ptrLibrary);
  // Subtract the underlying discrete vector of  the average function (uHat_B in paper)
  // from each column of W to create WHat
  LinearOperator<double> WHat = generateFluctuationMatrix(W);

  // Temporarily output the size of WHat
  cout << "WHat is " << WHat.range() << "x" << WHat.domain() << endl;

  // Create the POD object and obtain the POD basis functions
  POD_SVD pod(WHat, velocityDS, verbosity);
  pod.calculateSVD();
  pod.calculateBasisFunctions();
  
  // Define the relative path to where to store the results of the POD
  string POD_DataDir = "Results/tFinal"+Teuchos::toString(int(tFinal))+"sec/ReValues";
  // Format the displayed ReValues information system() can't create directories with
  // parenthesis in the name
  string ReValuesStr = "{";
  for(int ReIndex=0;ReIndex<ReValues.length()-1;ReIndex++)
      ReValuesStr = ReValuesStr + Teuchos::toString(ReValues[ReIndex]) + ",";
  ReValuesStr = ReValuesStr + Teuchos::toString(ReValues[ReValues.length()-1]) + "}";
  POD_DataDir += ReValuesStr + "/nx" + Teuchos::toString(nx)
    + "nt" + Teuchos::toString(nSteps) + "/tol";
  std::ostringstream tolFileValue;
  tolFileValue << std::setprecision(precision) << tol;
  POD_DataDir = POD_DataDir + tolFileValue.str() + "/";
  Array<Expr> phi = pod.get_basis_functions(tol,POD_DataDir);
  cout << "Returned from get_basis_functions(tol,POD_DataDir) " << std::endl;


  timer.stop();
  cout << "For tFinal = " << tFinal << " sec, Re = " << Re << ", nx = " << nx << ", nSteps = "
       << nSteps << ", and tol = " << tol << endl;
  cout << "runtime = " << timer.totalElapsedTime() << endl << endl << endl;

  Sundance::finalize(); 

}
