#include "Sundance.hpp"
#include "PlayaSerialVectorType.hpp" // For SerialVectorType
#include "PlayaDenseSerialMatrix.hpp" // For DenseSerialMatrix

//Viento
#include "VientoSnapshotIO.hpp"

double gradIP(Expr f, Expr g, Mesh mesh, QuadratureFamily quad)
{
	CellFilter interior = new MaximalCellFilter();
	// mesh.spatialDim() returns n for nD
        int dim = mesh.spatialDim();

	// Define our differential operators; note Derivative(x=0)
	Expr grad = gradient(dim);

	FunctionalEvaluator IP = FunctionalEvaluator(mesh, Integral(interior, colonProduct( outerProduct(grad,f),outerProduct(grad,g)), quad));
	return (IP.evaluate());
}

double tensorIP(Expr f, Expr g, Expr h, Mesh mesh, QuadratureFamily quad)
{
	CellFilter interior = new MaximalCellFilter();
	// mesh.spatialDim() returns n for nD
	// Define grad operator
	Expr grad = gradient(mesh.spatialDim());

//      The first one is what KL and I did on 1/27; then I realized the indices were off
//	FunctionalEvaluator IP = FunctionalEvaluator(mesh, Integral(interior, h*((f*grad)*g),quad));
	FunctionalEvaluator IP = FunctionalEvaluator(mesh, Integral(interior, f*((h*grad)*g),quad));
	return (IP.evaluate());
}

LinearOperator<double> makeA(Teuchos::Array<Expr> phi, Mesh mesh, QuadratureFamily quad, VectorType<double> vecType)
{
// Cause I don’t know how to get the type from a vector object
//	Playa::VectorType<double> vecType = phi[0].getVector().space().createMember();

	// Create an RxR matrix
	Playa::VectorSpace<double> domain = vecType.createEvenlyPartitionedSpace(Playa::MPIComm::world(), phi.length());
	Playa::VectorSpace<double> range = vecType.createEvenlyPartitionedSpace(Playa::MPIComm::world(), phi.length());
	Teuchos::RCP<Playa::MatrixFactory<double> > mf = vecType.createMatrixFactory(domain, range);
	Playa::LinearOperator<double> A = mf->createMatrix();

	// Access the matrix as a DenseSerialMatrix
	Teuchos::RCP<Playa::DenseSerialMatrix> APtr 
		        = Teuchos::rcp_dynamic_cast<Playa::DenseSerialMatrix>(A.ptr());
	/* Access the matrix as an EpetraMatrix
	Teuchos::RCP<Playa::EpetraMatrix> APtr 
		        = Teuchos::rcp_dynamic_cast<Playa::EpetraMatrix>(A.ptr());*/
	double* data = APtr->dataPtr();

	// a_{ij} = data[i + R*j]
	for(int i = 0; i<APtr->numRows(); i++)
	{
		for(int r = 0; r<APtr->numCols(); r++)
			data[i + APtr->numRows()*r] = gradIP(phi[i], phi[r],mesh,quad);
	}

	return A;
}



/********************************************************************************
*
*
********************************************************************************/
LinearOperator<double> makeT(Teuchos::Array<Expr> phi, Mesh mesh, QuadratureFamily quad, VectorType<double> vecType)
{
// Cause I don’t know how to get the type from a vector object
//	Playa::VectorType<double> vecType = phi[0].getVector().space().createMember();

	// Create an RxR matrix
	Playa::VectorSpace<double> domain = vecType.createEvenlyPartitionedSpace(Playa::MPIComm::world(), phi.length());
	Playa::VectorSpace<double> range = vecType.createEvenlyPartitionedSpace(Playa::MPIComm::world(), phi.length());
	Teuchos::RCP<Playa::MatrixFactory<double> > mf = vecType.createMatrixFactory(domain, range);
	Playa::LinearOperator<double> A = mf->createMatrix();

	// Access the matrix as a DenseSerialMatrix
	Teuchos::RCP<Playa::DenseSerialMatrix> APtr 
		        = Teuchos::rcp_dynamic_cast<Playa::DenseSerialMatrix>(A.ptr());
	/* Access the matrix as an EpetraMatrix
	Teuchos::RCP<Playa::EpetraMatrix> APtr 
		        = Teuchos::rcp_dynamic_cast<Playa::EpetraMatrix>(A.ptr());*/
	double* data = APtr->dataPtr();

	// a_{ij} = data[i + R*j]
	for(int i = 0; i<APtr->numRows(); i++)
	{
		for(int r = 0; r<APtr->numCols(); r++)
			data[i + APtr->numRows()*r] = gradIP(phi[i], phi[r],mesh,quad);
	}

	return A;
}

















int main(int argc, char** argv)
{

try{
	int nx = 32;
	Sundance::setOption("nx",nx, "Number of elements along each axis");

	Sundance::init(&argc, &argv);

	// Define our mesh type
	MeshType meshType = new Sundance::BasicSimplicialMeshType();
	// Define the dimensions of our domain
	double xmin = 0.0;
	double xmax = 1.0;
	double ymin = 0.0;
	double ymax = 1.0;
	// Build the 2D mesh
	Sundance::MeshSource mesher = new Sundance::PartitionedRectangleMesher(xmin,xmax,nx,1,ymin,ymax,nx,1,meshType,0);
	Sundance::Mesh mesh = mesher.getMesh();

	// Define our coordinates
	Expr x = new CoordExpr(0,"x");
	Expr y = new CoordExpr(1,"y");

	// Define our differential operators
	Expr dx = new Derivative(0);
	Expr dy = new Derivative(1);
	Expr grad = List(dx,dy);
/*
	Playa::VectorType<double> vecType = new Playa::EpetraVectorType();
	int num_vars = 2;	
	// mesh.spatialDim() returns n for nD
        int dim = mesh.spatialDim();
	Playa::VectorSpace<double> domain = vecType.createEvenlyPartitionedSpace(Playa::MPIComm::world(), num_vars);
	Playa::VectorSpace<double> range = vecType.createEvenlyPartitionedSpace(Playa::MPIComm::world(), dim);
	Teuchos::RCP<Playa::MatrixFactory<Expr> > mf = vecType.createMatrixFactory(domain, range);
	Playa::LinearOperator<Expr> jacobian = mf->createMatrix();
*/

	// Create our test functions
//	Expr f = List(x, y);
//	Expr g = List(x*x, x+y);
	Expr f = List(sin(x),exp(x));
	Expr g = List(x*x, cos(x+y));

	// Create the IP = \int_domain f*g
	CellFilter interior = new MaximalCellFilter();
	QuadratureFamily quad4 = new GaussianQuadrature(4);
//	FunctionalEvaluator IP = FunctionalEvaluator(mesh, Integral(interior, colonProduct( outerProduct(grad,f),outerProduct(grad,g)), quad4));
//	double value = IP.evaluate();

	double exactValue = -0.635104;
	std::cout << "The value for < (grad*f), (grad*g) > is: " << gradIP(f,g,mesh,quad4) << std::endl
		  << "The value should be " << exactValue << std::endl;
	/*std::cout << "f(x,y) is " << f << std::endl
		  << "grad*f is " << (grad*f) << std::endl
		  << "g(x,y) is " << g << std::endl
		  << "grad*g is " << (grad*g) << std::endl
		  << "first component of g is " << g[0] << std::endl
		  << "jacobian of f is " << outerProduct(f,grad) << std::endl
		  << "? of f is " << outerProduct(grad,f) << std::endl
		  << "Testing order " << x*dx << std::endl
                  << "Frobenius? is " << colonProduct(outerProduct(grad,f),outerProduct(grad,g)) << std::endl;
	*/	   

/*
Trying to figure out how Sundance handles the dot product; i.e. if i give it u*v, does it do u^t v?
	Expr test1 = List(2,3);
	Expr test2 = List(List(1,2),List(3,4));
	std::cout << "test1 dot test 2 " << test1*test2 << std::endl;
*/

	// Begin testing tensorIP
	Expr h1 = List(x*x, 3*y + 2*x);
//	Expr h1 = List(List(x*x), List(3*y+2*x));
//	Expr h1 = List(List(x*x, 3*y+2*x));
	Expr h2 = List(x*x*x*x + 10*y, y*y*y -1);
	Expr h3 = List(x*y+y*y, 24*y);
	Expr t = h3*((h1*grad)*h2);
	//std::cout << "Here is t: " << t << std::endl;
	exactValue = 101.708;
	std::cout << "The result of tensorIP(1,2,3) is " << tensorIP(h1,h2,h3,mesh,quad4) << std::endl
		  << "The value should be " << exactValue << std::endl;;


	// Testing out Teuchos::Array<Sundance::Expr>
	Teuchos::Array<Sundance::Expr> array = Teuchos::tuple(f,g);
	std::cout << "array is " << array << std::endl
		  << "size of array is " << array.size() << std::endl
		  << "length of array is " << array.length() << std::endl;

	Playa::VectorType<double> vecType = new Playa::SerialVectorType();
	Playa::LinearOperator<double> A = makeA(array, mesh, quad4, vecType);
	std::cout << "Here is A: " << std::endl <<	 A << std::endl;

	Teuchos::Array<Playa::LinearOperator<double> > tensorArray;








}

catch(exception& e)
{
	Sundance::handleException(e);
}

Sundance::finalize();
return 0;
}
