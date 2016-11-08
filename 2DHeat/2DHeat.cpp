/*
 * The purpose of this code is to solve the 2D Heat equation outlined by the jpeg in this folder.
 * This will give me experience with how to solve a problem in Sundance
 * Date 10/23/2016
 * Author: Simon  Rush
 *
*/

//#include "PlayaVectorType.hpp" // For VectorType, VectorSpace
#include "Sundance.hpp"

// Things necessary for mathematica generated fs to be understood
const double Pi = 4.0*atan(1.0);
Expr Cos(const Expr& x) {return cos(x);}
Expr Sin(const Expr& x) {return sin(x);}
Expr Power(const Expr& x, const double& p) {return pow(x,p);}
Expr Power(const double& x, const Expr& p) {return exp(p*log(x));}
const double E = exp(1.0);
Expr Sqrt(const Expr& x) {return sqrt(x);}
double Sqrt(const double& x) {return sqrt(x);}



int main (int argc, char *argv[])
{      
	/*
	 * This particular problem's error seems to be controlled by the number of timesteps
	 * nt = 100: .136
	 * nt = 400: .005
	 * nt = 800: .001
	*/
	int nx = 16;
	int verbosity = 0;
	int nt = 400;


      	//Use this options from the command line via ./executable --"name"=value
      	Sundance::setOption("nx",nx,"Number of elements along each axis");
      	Sundance::setOption("verbosity",verbosity,"Controls the level of detail the program prints to screen about what it is doing");
	Sundance::setOption("nt", nt, "Number of time steps used in solving the problem");

      	// Initialization steps (MPI, etc)
      	Sundance::init(&argc,&argv);

	// Define a coordiante
	Sundance::Expr x = new Sundance::CoordExpr(0);
	Sundance::Expr y = new Sundance::CoordExpr(1);
     	/* Represent the time variable as a parameter expression, NOT as
      	 * a double variable. The reason is that we need to be able to update
      	 * the time value without rebuilding expressions. 
	*/
     	Sundance::Expr t = new Sundance::Parameter(0.0);
	Sundance::Expr tPrev = new Sundance::Parameter(0.0);




	// Define our some parameters for our problem, the heat equation in 2D
	double xmin = -1.0;
	double xmax = 1.0;
	double ymin = -1.0;
	double ymax = 1.0;
	double tInit = 0.0;
	double tFinal = 8.0*Pi;
	double deltaT = (tFinal-tInit)/nt;
	int k = 3;
	double a = 0.75;
	double gamma = 2.0*Pi;
	Sundance::Expr Xp1 = a*((1.0+cos(t))/2.0)*cos(gamma*t);
	Sundance::Expr Xp2 = a*((1.0+cos(t))/2.0)*sin(gamma*t);

	
	// Define our vector type 
	Playa::VectorType<double> vecType = new Playa::EpetraVectorType();
	

      	// Build the mesh and domain
	Sundance::MeshType meshType = new Sundance::BasicSimplicialMeshType();
	// npx*npy = np (I am running serial, so 1)
      	Sundance::MeshSource mesher = new Sundance::PartitionedRectangleMesher(xmin, xmax, nx,1,ymin,ymax,nx,1, meshType,verbosity);
      	Sundance::Mesh mesh = mesher.getMesh();
	
	// Create a variable to store the number of dimensions of the problem
	// mesh.spatialDim() returns n for nD
        int dim = mesh.spatialDim();

	// Set up differential operators; expects n for nD
	Sundance::Expr grad = gradient(dim);

	// Filter subtype MaximalCellFilter selects all cells having dimension equal to the spatial
	// dimension of the mesh. 
      	Sundance::CellFilter interior = new Sundance::MaximalCellFilter();
    	// DimensionalCellFilter expects dimesnions to start at 0. It returns 
      	// the boundaries of the mesh; i.e. the edges of a square. nD is ( n-1)
      	Sundance::CellFilter sides = new Sundance::DimensionalCellFilter(dim-1);
      	Sundance::CellFilter north = sides.subset(new Sundance::CoordinateValueCellPredicate(1,ymax));
      	Sundance::CellFilter east = sides.subset(new Sundance::CoordinateValueCellPredicate(0,xmax));
      	Sundance::CellFilter south = sides.subset(new Sundance::CoordinateValueCellPredicate(1,ymin));
      	Sundance::CellFilter west = sides.subset(new Sundance::CoordinateValueCellPredicate(0,xmin));

      	// BasisFamily used to express our solutions
      	Sundance::BasisFamily basis = new Sundance::Lagrange(2); // 2nd order PWQL (Piece-Wise Quadratic Lagrange)

	// Test Functions
	Sundance::Expr v = new Sundance::TestFunction(basis, "v");
/*
	int numTest = 1;
	Teuchos::Array<Expr> v(numTest);
	for(int i=0; i<v.size(); i++)
		v[i] = new TestFunction(basis, "v[" + Teuchos::toString(i) + "]");

	Sundance::Expr vlist = new Sundance::ListExpr(v);
*/

	// Unknown Functions
 	Sundance::Expr u = new Sundance::UnknownFunction(basis, "u");
/*
	int numUnknown = 1;
	Teuchos::Array<Expr> u(numUnknown);
	for(int i=0; i<u.size(); i++)
		u[i] = new UnknownFunction(basis, "u[" + Teuchos::toString(i) + "]");
	Sundance::Expr ulist = new Sundance::ListExpr(u);
*/


      // Evaluation scheme. The parameter says for what degree of polynomials
      // it will be exact for
      Sundance::QuadratureFamily quad = new Sundance::GaussianQuadrature(4);


	Sundance::DiscreteSpace discSpace(mesh,basis,vecType);
	Sundance::Expr uExact = exp(-k*( (x-Xp1)*(x-Xp1) + (y-Xp2)*(y-Xp2) ))*(1-x*x)*(1-y*y);
	Sundance::Expr q = (Power(E,-3*Power(x - 
          (3*(1 + Cos(t))*Cos(2*Pi*t))/8.,2) - 
       3*Power(y - (3*(1 + Cos(t))*Sin(2*Pi*t))/
           8.,2))*(64*
        (8 - 37*Power(x,2) + 18*Power(x,4) + 
          (-37 + 66*Power(x,2) - 18*Power(x,4))*
           Power(y,2) - 
          18*(-1 + Power(x,2))*Power(y,4)) - 
       648*(-1 + Power(x,2))*(-1 + Power(y,2))*
        Power(Cos(t/2.),4) + 
       9*(-1 + Power(x,2))*(-1 + Power(y,2))*
        Sin(t)*(3 + 3*Cos(t) - 8*x*Cos(2*Pi*t) - 
          8*y*Sin(2*Pi*t)) + 
       288*Power(Cos(t/2.),2)*
        ((-1 + Power(y,2))*
           (-(Pi*y) + 
             x*(-8 + 6*Power(x,2) + Pi*x*y))*
           Cos(2*Pi*t) - 
          (-1 + Power(x,2))*
           (8*y - 6*Power(y,3) + 
             Pi*x*(-1 + Power(y,2)))*Sin(2*Pi*t)))
     )/32.;



	Sundance::Expr qPrev = (Power(E,-3*Power(x - 
          (3*(1 + Cos(tPrev))*Cos(2*Pi*tPrev))/8.,
         2) - 3*Power(y - 
          (3*(1 + Cos(tPrev))*Sin(2*Pi*tPrev))/8.,
         2))*(64*(8 - 37*Power(x,2) + 
          18*Power(x,4) + 
          (-37 + 66*Power(x,2) - 18*Power(x,4))*
           Power(y,2) - 
          18*(-1 + Power(x,2))*Power(y,4)) - 
       648*(-1 + Power(x,2))*(-1 + Power(y,2))*
        Power(Cos(tPrev/2.),4) + 
       9*(-1 + Power(x,2))*(-1 + Power(y,2))*
        Sin(tPrev)*(3 + 3*Cos(tPrev) - 
          8*x*Cos(2*Pi*tPrev) - 
          8*y*Sin(2*Pi*tPrev)) + 
       288*Power(Cos(tPrev/2.),2)*
        ((-1 + Power(y,2))*
           (-(Pi*y) + 
             x*(-8 + 6*Power(x,2) + Pi*x*y))*
           Cos(2*Pi*tPrev) - 
          (-1 + Power(x,2))*
           (8*y - 6*Power(y,3) + 
             Pi*x*(-1 + Power(y,2)))*
           Sin(2*Pi*tPrev))))/32.;


	//proj() expects the discrete space to have an EpetraVectorType
	Sundance::L2Projector proj_u(discSpace, uExact);
	Sundance::Expr uPrev = proj_u.project();

	// Define what you want to integrate
	Sundance::Expr integrand = ( (u-uPrev)/deltaT )*v + 0.5*(grad*u + grad*uPrev)*(grad*v) - 0.5*(q + qPrev)*v;
	Sundance::Expr eqn = Integral(interior,integrand,quad);

	// Define Dirichlet BC
	Sundance::Expr h = new Sundance::CellDiameterExpr();
	Sundance::Expr bc = Sundance::EssentialBC(north+east+south+west, v*u,quad);

	std::cout << "Defining Problem........" << std::endl;
	Sundance::LinearProblem prob(mesh,eqn,bc,v,u,vecType);


	Sundance::LinearSolver<double> solver = Sundance::LinearSolverBuilder::createSolver("amesos.xml");

	// Set up the directory to write the results
	const int sys_error = system( ("rm -fr results/nx"+Teuchos::toString(nx)+"/nt"+std::to_string(nt) ).c_str() );
	if(-1==sys_error)
	{
		std::cout << "Error removing results/nx" << nx << "/nt" << nt << std::endl;
		exit(1);
	}


	string outputDir = "results/nx"+Teuchos::toString(nx);
	const int nx_dir_error = system( ("mkdir -p " + outputDir).c_str() );	
	if (-1 == nx_dir_error)
	{
  		  printf("Error creating nx directory!\n");
  	 	  exit(1);
	}
/*
	const int nx_dir_error = mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if (-1 == nx_dir_error)
	{
  		  printf("Error creating nx directory!\n");
  	 	  exit(1);
	}
*/	

	outputDir = outputDir + "/nt"+std::to_string(nt);
	const int nt_dir_error = system( ("mkdir " + outputDir).c_str() );
	if (-1 == nt_dir_error)
	{
  		  printf("Error creating nt directory!\n");
  	 	  exit(1);
	}
/*
	const int nt_dir_error = mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if (-1 == nt_dir_error)
	{
  		  printf("Error creating nt directory!\n");
  	 	  exit(1);
	}	
*/
	

	Sundance::FieldWriter writer0 = new Sundance::VTKWriter(outputDir + "/2DHeat-0");
	writer0.addMesh(mesh);
	writer0.addField("T", new ExprFieldWrapper(uPrev[0]));
	writer0.write();

	double error = 0.0;
	
	for(int i = 0; i < (tFinal-tInit)/deltaT; i++)
	{
		t.setParameterValue((i+1)*deltaT);
		tPrev.setParameterValue(i*deltaT);
		std::cout << "t=" << (i+1)*deltaT << std::endl;
		Sundance::Expr uNext = prob.solve(solver);
		error = Sundance::L2Norm(mesh, interior, uExact-uNext, quad);
		std::ostringstream oss;
		oss << outputDir + "/2DHeat-" << i+1;
		Sundance::FieldWriter writer = new Sundance::VTKWriter(oss.str());
		writer.addMesh(mesh);
		writer.addField("T", new Sundance::ExprFieldWrapper(uNext[0]));
		writer.write();

		updateDiscreteFunction(uNext, uPrev);
	}


		error = Sundance::L2Norm(mesh, interior, uExact-uPrev, quad);
		std::cout << "error norm=" << error << std::endl;
		/* 
		 * Only for error checking. If the functional norm is large
		 * then the error norm doesn't have to be tiny
		 * double funcNorm = Sundance::L2Norm(mesh, interior, uExact, quad);
		 * std::cout << "Function norm=" << funcNorm << std::endl;
		*/
	


	Sundance::finalize();
	return 0;
}
