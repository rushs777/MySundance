#include "Sundance.hpp"

// Standard Library Functions
using std::vector;
using std::cout;
using std::endl;

// Namespaces
using namespace Teuchos;
using namespace Playa;

int main(int argc, char *argv[])
{
  try
    {
      Time timer("total");
      timer.start();

      // Sundance Options
      int nSteps = 100;
      Sundance::setOption("nSteps",nSteps,"Number of time steps taken");

      // Initialization steps (MPI, etc)
      Sundance::init(&argc,&argv);

      // Construct the framework for the optimization
      double tInit = 0.0;
      double T = 1.0;
      double deltaT = T/((double) nSteps);
      double pmin = -5.0;
      double pmax = 5.0;
      // Represent the time variable as a Sundance::Parameter so that the time value
      // can be updated and the problem solved without having to rebuild the problem
      Expr t = new Sundance::Parameter(tInit);
      Expr tPrev = new Sundance::Parameter(tInit);

      double alpha = 6.0;
      Expr xExact = alpha*exp(t);

      // Build the time mesh
      MeshType timeMeshType = new BasicSimplicialMeshType();
      MeshSource timeMesher = new PartitionedLineMesher(tInit, T, nSteps, timeMeshType, 0);
      Mesh timeMesh = timeMesher.getMesh();

      // Filter subtype MaximalCellFilter selects all cells having dimension equal to the
      // spatial dimension of the mesh
      CellFilter interior = new MaximalCellFilter();
      // DimensionalCellFilter 0=vertices, 1=edges, 2=triangles, 3=tetrahedra
      CellFilter points = new DimensionalCellFilter(0);
      CellFilter tStart = points.coordSubset(0, tInit);
      CellFilter tEnd = points.coordSubset(0, T);

      // Define our DiscreteSpace(mesh, basis, vector type)
      VectorType<double> vecType = new EpetraVectorType();
      BasisFamily xbasis = new Lagrange(2);
      DiscreteSpace timeDS(timeMesh, xbasis, vecType);

      // Define our Unknown and Test functions
      Expr v = new TestFunction(xbasis,"v");
      Expr x = new UnknownFunction(xbasis,"x");

      // Define the quadrature rule to use for solving the problem
      QuadratureFamily quad = new GaussianQuadrature(4);

      // Create xPrev through projection
      L2Projector projX(timeDS, xExact);
      Expr xPrev = projX.project();

      // Set up the LinearProblem
      Expr integrand = (1.0/deltaT)*(x-xPrev)*v - 0.5*(x+xPrev)*v;
      Expr eqn = Integral(interior,integrand,quad);
      Expr BC = EssentialBC(tStart,alpha,quad); // x(tStart) = alpha

      cout << "Defining Problem......." << endl;
      LinearProblem prob(timeMesh,eqn,BC,v,x,vecType);

      cout << "Solving Problem........" << endl;
      LinearSolver<double> solver = LinearSolverBuilder::createSolver("amesos.xml");
      Array<double> error(nSteps);
      Array<Expr> xApprox(nSteps+1); // Holds (t,x) parings as xApprox[t]
      xApprox[0] = alpha;
      for(int i=0; i < nSteps; i++)
	{
	  t.setParameterValue( (i+1)*deltaT ); //time value to solve for
	  tPrev.setParameterValue( i*deltaT );
	  cout << "For t=" << t.getParameterValue();
	  xApprox[i+1] = prob.solve(solver);
	  error[i+1] = L2Norm(timeMesh, interior, xExact - xApprox[i+1], quad);
	  cout << " the error is " << error[i+1] << endl;
	  cout << "xExact: " << L2Norm(timeMesh, interior, xExact, quad) << "\t xApprox: " << L2Norm(timeMesh, interior, xApprox[i+1], quad) << endl;

	  updateDiscreteFunction(xApprox[i+1], xPrev);
	}
      
      
      
      

      Sundance::finalize();
    }
  catch(std::exception& e)
    {
      Out::root() << "Caught exception: " << e.what() << endl;
      return -1;
    }
  


}
  
