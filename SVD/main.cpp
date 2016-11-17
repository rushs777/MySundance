#include "Teuchos_GlobalMPISession.hpp" //Handles initializing and finalizing MPI. Also defines rcp_dynamic_cast
#include "PlayaVectorType.hpp" // For VectorType, VectorSpace
#include "PlayaSerialVectorType.hpp" // For SerialVectorType
//#include "/home/sirush/PhDResearch/DenseMatrixMatrixProduct/PlayaDenseMatrixMatrixProduct.hpp"
#include "PlayaDenseSerialMatrix.hpp" // For DenseSerialMatrix
#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaSimpleTransposedOpDecl.hpp"
#include "VientoSnapshotIO.hpp"

#include "PlayaSVD.hpp"
// For SVD's main.cpp
#include "Sundance.hpp"

/*
#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp" // Needed to create an empty LinearOperator object
#include "PlayaSimpleTransposedOpImpl.hpp" //Needed if I want to take transposes

#include "PlayaSimpleComposedOpImpl.hpp" //Needed for matrix, matrix multiply via *
#endif
*/

// Is this how I want to do this?
const double PI = 4.0*atan(1.0);

// Temporary addition until I can get Viento to work right
/*
#include <fstream>
#include <iomanip>
using std::ifstream;
using std::to_string;
using std::setw;
using std::setprecision;
Expr readSnap(const string& filename, int tag, const Mesh& mesh);
LinearOperator<double> snapshotToMatrix(const string& filename, int maxTagNum, const Mesh& mesh);
*/
int main (int argc, char *argv[])
{      
	int nx = 1;
	int verbosity = 0;
	int PhiCol = 0;

      //Use this options from the command line via ./executable --"name"=value
      Sundance::setOption("nx",nx,"Number of elements along each axis");
      Sundance::setOption("verbosity",verbosity,"Controls the level of detail the program prints to screen about what it is doing");
	Sundance::setOption("PhiCol",PhiCol,"Controls what column of Phi is printed out");
      // Initialization steps (MPI, etc)
      Sundance::init(&argc,&argv);

	// Define our vector type 
	Playa::VectorType<double> vecType = new Playa::SerialVectorType();
	// Define our mesh type
	Sundance::MeshType meshType = new Sundance::BasicSimplicialMeshType();

// Starting Here and going forward, all of this code can be removed by the fact that POD now creates
// the mass matrix from the mesh only
/*
	// Define our some parameters for our problem, the heat equation in 2D
	double xmin = -1.0;
	double xmax = 1.0;
	double ymin = -1.0;
	double ymax = 1.0;	

      	// Build the mesh and domain
	// npx*npy = np (I am running serial, so 1)
      	Sundance::MeshSource mesher = new Sundance::PartitionedRectangleMesher(xmin, xmax, nx,1,ymin,ymax,nx,1, meshType,verbosity);
      	Sundance::Mesh mesh = mesher.getMesh();
	
	// Create a variable to store the number of dimensions of the problem
	// mesh.spatialDim() returns n for nD
        int dim = mesh.spatialDim();
	std::cout << "dim is " << dim << std::endl;

	// Filter subtype MaximalCellFilter selects all cells having dimension equal to the spatial
	// dimension of the mesh. 
      	Sundance::CellFilter interior = new Sundance::MaximalCellFilter();


      // BasisFamily used to express our solutions
      Sundance::BasisFamily basis = new Sundance::Lagrange(2); // 2nd order PWQL (Piece-Wise Quadratic Lagrange)

	// Test Functions
	int numTest = 1;
	Teuchos::Array<Expr> v(numTest);
	for(int i=0; i<v.size(); i++)
		v[i] = new TestFunction(basis, "v[" + Teuchos::toString(i) + "]");

	
	Sundance::Expr vlist = new Sundance::ListExpr(v);

	// Unknown Functions
	int numUnknown = 1;
	Teuchos::Array<Expr> u(numUnknown);
	for(int i=0; i<u.size(); i++)
		u[i] = new UnknownFunction(basis, "u[" + Teuchos::toString(i) + "]");

	
	Sundance::Expr ulist = new Sundance::ListExpr(u);

      // Evaluation scheme. The parameter says for what degree of polynomials
      // it will be exact for
      Sundance::QuadratureFamily quad = new Sundance::GaussianQuadrature(4);

	// Define what you want to integrate
	Sundance::Expr integrand = vlist*ulist;
	Sundance::Expr eqn = Integral(interior,integrand,quad);
	
	// Define Empty BC
	Sundance::Expr bc;// since I want the mass matrix

	std::cout << "Defining Problem........" << std::endl;
	Sundance::LinearProblem prob(mesh,eqn,bc,vlist,ulist,vecType);

	std::cout << "Getting S..........." << std::endl;
	Playa::LinearOperator<double> S = prob.getOperator();
	std::cout << "S is a " << S.range().dim() << " by " << S.domain().dim() << std::endl;
	//std::cout << "Here is S" << std::endl << S << std::endl;

	std::cout << "Writing mesh........." << std::endl;
	Sundance::FieldWriter meshWriter = new Sundance::ExodusWriter("testmesh"); //give it filename, not filename.ext 
	meshWriter.addMesh(mesh);
	meshWriter.write(); //adds the .exo extension

	//	std::cout << "nPts = " << mesh.numCells(0) << std::endl
	//		  << "nodePosition(0) = " << mesh.nodePosition(0) << std::endl; 
*/
// End of the code to get S/mesh, which is no longer necessary

	std::cout << "Reading mesh.........." << std::endl;
	Sundance::MeshSource meshReader = new Sundance::ExodusMeshReader("cyl-0", meshType, verbosity);
	Sundance::Mesh mesh = meshReader.getMesh();

/*	int dim = 0;
	std::cout << "Here is the number of cells along the x-axis: " << mesh.numCells(dim) << std::endl
		  << "Here is the number of cells along the y-axis: " << mesh.numCells(1) << std::endl;

*/
	// Set up W for cyl-0.exo
	Playa::LinearOperator<double> W = snapshotToMatrix("results/st-p", 3, mesh);

	// Get the POD
	Playa::LinearOperator<double> Alpha;
	Playa::LinearOperator<double> Phi;
	Playa::Vector<double> lambda;

	// W and mesh need to be defined
	POD(W,lambda,Alpha,Phi,mesh);
	std::cout << "POD finished" << std::endl;


	//Start looking at svd.cpp

	// Set up the snapshot matrix W
/*	int gridPts = 9; // rows
	int timeSteps = 10; // cols
	Playa::VectorSpace<double> domain = vecType.createEvenlyPartitionedSpace(Playa::MPIComm::world(),timeSteps);
	Playa::VectorSpace<double> range = vecType.createEvenlyPartitionedSpace(Playa::MPIComm::world(),gridPts);
	Teuchos::RCP<Playa::MatrixFactory<double> > mfW = vecType.createMatrixFactory(domain, range);
	Playa::LinearOperator<double> W = mfW->createMatrix();
	Teuchos::RCP<Playa::DenseSerialMatrix> PtrW = Teuchos::rcp_dynamic_cast<Playa::DenseSerialMatrix>(W.ptr());
	// Initialize W
	for(int row = 0; row<PtrW->numRows(); row++)
		for(int col = 0; col<PtrW->numCols(); col++)
			PtrW->dataPtr()[row+PtrW->numRows()*col] = col+PtrW->numCols()*row;
*/
	/* Checking behavior of W*e_j
	Playa::Vector<double> ej = domain.createMember();
	ej.zero();
	ej[0]=1.0;
	Playa::Vector<double> result;
	W.apply(ej,result);
	std::cout << "Here is W*e_j " << std::endl << result << std::endl;*/



	// Create lambda, Alpha, and Phi
	//Playa::LinearOperator<double> Alpha;
	//Playa::LinearOperator<double> Phi;
	//Playa::Vector<double> lambda;

	// W and mesh need to be defined
	//POD(W,lambda,Alpha,Phi,mesh);

	//std::cout << "Here is S: " << std::endl << S << std::endl;
	//std::cout << "Here is W: " << std::endl << W << std::endl;
	//std::cout << "Here is lambda " << std::endl << lambda << std::endl;
	//std::cout << "Here is Alpha " << std::endl << Alpha << std::endl;
	//std::cout << "Here is Phi " << std::endl << Phi << std::endl;
	
	// Test the values returned from POD
	//std::cout << "Here is lambda " << std::endl << lambda << std::endl;
	//std::cout << "Here is Alpha " << std::endl << Alpha << std::endl;
	//std::cout << "Here is Phi " << std::endl << Phi << std::endl;
	//Teuchos::RCP<Playa::DenseSerialMatrix> PhiPtr = 						
		Teuchos::rcp_dynamic_cast<Playa::DenseSerialMatrix>(Phi.ptr());

	//std::cout << "Here is Phi[[All," << PhiCol << "]]" << std::endl;
	//for(int row = 0; row < PhiPtr->numRows(); row++)
		//std::cout << PhiPtr->dataPtr()[row+PhiPtr->numRows()*PhiCol] << std::endl;

/*
test code for snapshotToMatrix

string filename = "results/st-v";
MeshSource mesher = new ExodusMeshReader("cyl-0", meshType);
int tagSize = 5;

	Playa::LinearOperator<double> B = snapshotToMatrix(filename, tagSize, mesher.getMesh() );

	std::cout << "Here's B: " << std::endl << B << std::endl;
*/


/*
	Vector<double> error = getDiscreteFunctionVector(ptest);
	Playa::Vector<double> ej = d.createMember();
	ej.zero();
	ej[0]=1.0;
	Playa::Vector<double> result;
	A.apply(ej,result);
	std::cout << "Size of ptest: " << error.dim() << std::endl
                  << "Size of result: " << result.dim() << std::endl;

	std::cout << "Attempting subtract the two vectors " << std::endl;
	error -= result;

	Playa::Vector<double> ej = d.createMember();
	ej.zero();
	ej[0]=1.0;
	Playa::Vector<double> result;
	A.apply(ej,result);
	error -= result;
	std::cout << "Double checking that this puts " << filename + to_string(tag) + ".vec into the matrix column: " << error.max() << std::endl;
*/
	
	return 0;
}


