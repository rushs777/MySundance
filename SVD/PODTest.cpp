#include "Sundance.hpp" // Contains a variety of header files. See website for full list

//Needed for dense matrices
#include "PlayaSerialVectorType.hpp"
#include "PlayaDenseSerialMatrix.hpp"

//Viento
#include "VientoSnapshotIO.hpp"

//Local Files
#include "PlayaSVD.hpp"

using std::cout;
using std::endl;

int main (int argc, char *argv[])
{
	// Define Sundance Options
	int nx = 16;
	int debug = 1;

	Sundance::setOption("nx",nx,"Number of elements along each axis");
	Sundance::setOption("debugLevel",debug,"Controls the level of debugging information that is printed out");

	// Initialize Sundance
	Sundance::init(&argc, &argv);

	// Define our vector type
	//Playa::VectorType<double> vecType = new Playa::SerialVectorType();
	Playa::VectorType<double> vecType = new Playa::EpetraVectorType();

	// Define our mesh type
	Sundance::MeshType meshType = new Sundance::BasicSimplicialMeshType();
	// Define the dimensions of our domain
	double xmin = -1.0;
	double xmax = 1.0;
	double ymin = -1.0;
	double ymax = 1.0;
	double zmin = -1.0;
	double zmax = 1.0;
	// Build the 2D mesh
	// npx*npy = np; serial means these are all 1
	Sundance::MeshSource mesher2D = new Sundance::PartitionedRectangleMesher(xmin,xmax,nx,1,ymin,ymax,nx,1,meshType,debug-1);
	Sundance::Mesh mesh2D = mesher2D.getMesh();
	// Extract a 3D mesh from the 2D mesh
	Sundance::ExtrusionMeshTransformation extrusion(zmin,zmax,nx,meshType);
	Sundance::Mesh mesh3D = extrusion.apply(mesh2D);

	if(debug >=2)
	{
		SUNDANCE_ROOT_MSG1(debug, "Writing mesh......");
		Sundance::FieldWriter meshWriter3D = new Sundance::ExodusWriter("cube3D");
		meshWriter3D.addMesh(mesh3D);
		meshWriter3D.write(); // adds the .exo extension
		SUNDANCE_ROOT_MSG2(debug, "Mesh written........");
	}

	// Define our coordinates
	Expr x = new CoordExpr(0, "x");
	Expr y = new CoordExpr(1, "y");
	Expr z = new CoordExpr(2, "z");

	int n = 3;
	Teuchos::Array<Expr> psi(n);
	Teuchos::Array<Expr> phi(n);
	for(int count=0; count<psi.size(); count++)
	{
		psi[count] = sin( (double)count*x)*cos( (double)count*y)*z*z;
		Expr psi_k = List(0.0, 0.0, psi[count]);
		phi[count] = curl(psi_k);
	}

	int numVertices = mesh3D.numCells(0);
	int numEdges = mesh3D.numCells(1);
	int numNodes = numVertices + numEdges;
	cout << "Number of vertices: " << numVertices << endl
	     << "Number of edges: " << numEdges << endl
	     << "Number of nodes: " << numNodes << endl;


	// BasisFamily used to express our solutions
	Array<Sundance::BasisFamily> basis = List(new Sundance::Lagrange(2), new Sundance::Lagrange(2), new Sundance::Lagrange(2)); // 2nd order Piece-Wise Quad Lagrange
	//Sundance::BasisFamily basis = new Sundance::Lagrange(2); // 2nd order Piece-Wise Quad Lagrange

	// Create our time parameter
	double tInit = 0.0;
	Expr t(new Sundance::Parameter(tInit));
	double tFinal = 1.0;
	int Nt = 2;
	double deltaT = (tFinal-tInit)/Nt;
	

	// Create the function u
	Expr u = List(0.0, 0.0, 0.0);
	Teuchos::Array<double> c_n(n);
	Teuchos::Array<double> w_n(n);
	Teuchos::Array<double> d_n(n);
/*
	VectorSpace<double> vecSpace = vecType.createSpace(n);
	Vector<double> c_n = vecSpace.createMember();
	Vector<double> w_n = vecSpace.createMember();
	Vector<double> d_n = vecSpace.createMember();
*/
	for(int count = 0.0; count < n; count++)
	{
		c_n[count] = 1.0/(count+1);
		w_n[count] = 2.0/(3*count+1);
		d_n[count] = 4.0/(count*count+1);
		u = u + phi[count]*c_n[count]*cos(w_n[count]*t+d_n[count]);
	}

	// Build our L2Projector for evaluating a function on the mesh
	Sundance::DiscreteSpace ds(mesh3D, basis, vecType);
	Expr evaluated_u;

	SUNDANCE_ROOT_MSG1(debug, "Projecting u onto the mesh");
	for(int count=0; count<=Nt; count++)
	{
		SUNDANCE_ROOT_MSG1(debug,"time step "+std::to_string(count)+" of "+std::to_string(Nt));
		L2Projector projector(ds,u);
		evaluated_u = projector.project();
		writeSnap("results/cube3DTest",count,evaluated_u);
		t.setParameterValue(t.getParameterValue()+deltaT);
	}

	// Read the snapshots into a matrix
	Playa::LinearOperator<double> W = snapshotToMatrix("results/cube3DTest", Nt, mesh3D);
	cout << "Size of W: " << W.range().dim() << " by " << W.domain().dim() << endl;

	// Perform the POD of the matrix W
	Playa::LinearOperator<double> Alpha;
	Playa::LinearOperator<double> Phi;
	Playa::Vector<double> lambda;

	// W and mesh need to be defined
	SUNDANCE_ROOT_MSG1(debug, "Entering POD");
	POD(W,lambda,Alpha,Phi,mesh3D, debug);
	SUNDANCE_ROOT_MSG2(debug, "POD finished");

/*
	Teuchos::Array<Expr> evaluated_phi(n);

	for(int count=0; count<evaluated_phi.size(); count++)
	{
		L2Projector phiProjector(ds,phi[count]);
		evaluated_phi[count] = phiProjector.project();
		writeSnap("results/cube3DTest",count,evaluated_phi[count]);
	}

	Expr testfn = x + y +z;
	cout << "Basis size: " << basis.dim() << endl
 	     << "Expr size: " << testfn.size() << endl;
	L2Projector testfnProj(ds, testfn);
	Expr evaluatedfn = testfnProj.project();
	//setDiscreteFunctionVector(testfn, getDiscreteFunctionVector(evaluatedfn));
	writeSnap("cube3DTest",99,evaluatedfn);
*/
	


	


	return 0;
}
