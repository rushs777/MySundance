#include "Sundance.hpp" // Contains a variety of header files. See website for full list
#include "PlayaSerialEpetraAdapter.hpp" // Allows me to convert between Epetra and Serial Vectors

//Needed for dense matrices
#include "PlayaSerialVectorType.hpp"
#include "PlayaDenseSerialMatrix.hpp"

//Viento
#include "VientoSnapshotIO.hpp"

//Local Files
#include "PlayaSVD.hpp"

using std::cout;
using std::endl;

CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-2.0*M_PI) < 1.0e-10;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-2.0*M_PI) < 1.0e-10;})

CELL_PREDICATE(CornerPointTest, {return fabs(x[1]) < 1.0e-10 && fabs(x[0])<1.0e-10;})

int main (int argc, char *argv[])
{
	// Define Sundance Options
	int nx = 32;
	int verbosity = 1;
	int redoMesh = 0;
	double tFinal = 1.0;
	int Nt = 32;

	Sundance::setOption("nx",nx,"Number of elements along each axis");
	Sundance::setOption("verbosity",verbosity,"Controls the level of debugging information that is printed out");
	Sundance::setOption("redoMesh",redoMesh,"If the mesh needs to be rewritten, change to 1");
	Sundance::setOption("tFinal",tFinal,"The final time value");
	Sundance::setOption("Nt",Nt,"The number of timesteps to go from tInit to tFinal");

	// Initialize Sundance
	Sundance::init(&argc, &argv);

	// Define our vector type
	Playa::VectorType<double> vecType = new Playa::EpetraVectorType();

	// Define our mesh type
	Sundance::MeshType meshType = new Sundance::BasicSimplicialMeshType();
	// Define the dimensions of our domain
	double xmin = 0.0;
	double xmax = 2.0*M_PI;
	double ymin = 0.0;
	double ymax = 2.0*M_PI;
	// Build the 2D mesh
	// npx*npy = np; serial means these are all 1
	Sundance::MeshSource mesher = new Sundance::PartitionedRectangleMesher(xmin,xmax,nx,1,ymin,ymax,nx,1,meshType,verbosity-1);
	Sundance::Mesh mesh = mesher.getMesh();

	if(redoMesh!=0)
	{
		Sundance::FieldWriter meshWriter = new Sundance::ExodusWriter("0x2pi_2D");
		meshWriter.addMesh(mesh);
		meshWriter.write(); // adds the .exo extension
	}

	// Define our coordinates
	Expr x = new CoordExpr(0,"x");
	Expr y = new CoordExpr(1,"y");

	// Create our time parameter
	double tInit = 0.0;
	Expr t(new Sundance::Parameter(tInit));
	double deltaT = (tFinal-tInit)/Nt;

	// Get some information about our mesh
	int numVertices = mesh.numCells(0);
	int numEdges = mesh.numCells(1);
	int numNodes = numVertices + numEdges;
	cout << "Number of nodes: " << numNodes << endl;

	// Define our answer
	Expr uEx = List(sin(x)*cos(y), -cos(x)*sin(y))*exp(-2*t);

	// Create a BasisFamily to express our solution
	Array<Sundance::BasisFamily> basis = List(new Sundance::Lagrange(2), new Sundance::Lagrange(2)); // 2nd order Piece-Wise Quad Lagrange
/*
	// Create a discrete space to project onto
	Sundance::DiscreteSpace ds(mesh, basis, vecType);
	Expr evaluated_uEx;

	SUNDANCE_ROOT_MSG1(verbosity, "Projecting uEx onto the mesh");
	for(int count=0;count<=Nt; count++)
	{
		SUNDANCE_ROOT_MSG1(verbosity,"time step "+std::to_string(count)+" of "+std::to_string(Nt)+"; t is "+std::to_string(t.getParameterValue()));
		L2Projector projector(ds,uEx);
		evaluated_uEx = projector.project();
		writeSnap("results/Noodle/-1x1_2DTest",count,evaluated_uEx);
		t.setParameterValue(t.getParameterValue()+deltaT);		
	}
*/

	// Read the snapshots into a matrix
	Playa::LinearOperator<double> W = snapshotToMatrix("Results/nx32/nt32/st-v", Nt, mesh);
	SUNDANCE_ROOT_MSG2(verbosity, "Size of W: " << W.range().dim() << " by " << W.domain().dim());

	// Perform the POD of the matrix W
	Playa::LinearOperator<double> Omega;
	Playa::LinearOperator<double> Phi;
	Playa::Vector<double> lambda;

	// W and mesh need to be defined
	SUNDANCE_ROOT_MSG1(verbosity, "Entering POD");
	POD(W,lambda,Omega,Phi,mesh,verbosity);
	SUNDANCE_ROOT_MSG2(verbosity, "POD finished");
	
	//Needed for the integral
	CellFilter interior = new MaximalCellFilter();
	QuadratureFamily quad4 = new GaussianQuadrature(4);
	//Wanted these for solving the NSE
	CellFilter bdry = new BoundaryCellFilter();
	CellFilter top = bdry.subset( new TopPointTest() );
	CellFilter right = bdry.subset( new RightPointTest() );
	CellFilter bottom = bdry.subset( new BottomPointTest() );
	CellFilter left = bdry.subset( new LeftPointTest() );
	CellFilter nodes = new DimensionalCellFilter(0);
	CellFilter corner = nodes.subset(new CornerPointTest());


	Sundance::DiscreteSpace ds(mesh, basis, vecType);
	Expr phi_r;

	// Looking at PlayaSVD.cpp, Phi is a DenseSerialMatrix
	Playa::Vector<double> ej = Phi.domain().createMember();
	Playa::Vector<double> alpha = serialToEpetra(Phi.domain().createMember());
	Playa::Vector<double> phi;
	Expr testfn = x+y;
	alpha.zero();


	
	for(int r = 0; r<=Nt; r++)
	{
		ej.zero();
		ej[r] = 1.0;
	//std::cout << "After zeroing out ej" << std::endl;
		Phi.apply(ej,phi);
	//std::cout << "After getting phi_r" << std::endl;
		
		phi_r = new DiscreteFunction(ds, serialToEpetra(phi));	
		FunctionalEvaluator L2NormEvaluator = FunctionalEvaluator(mesh, Integral(interior, uEx*phi_r, quad4));
		alpha[r] = L2NormEvaluator.evaluate();
		t.setParameterValue(t.getParameterValue()+deltaT);

		//Expr temp = Integral(interior, uEx*phi_r, quad4);
		//double value = evaluateIntegral(mesh, temp);		
	}

	std::cout << "Here is the integral: " << alpha << endl;



	











	return 0;
}
