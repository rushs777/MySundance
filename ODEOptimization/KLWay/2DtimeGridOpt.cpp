#include "Sundance.hpp"
#include "PlayaAmesosSolver.hpp"

// Local Files
#include "MathematicaConverter.hpp"


// Standard Library Functions
using std::cout;
using std::endl;
using std::vector;

int main(int argc, char *argv[]) 
{
  try
    {
      int verb = 0;      
      Sundance::setOption("verbosity", verb, "verbosity level");
      
      int nx = 128;
      Sundance::setOption("nx", nx, "Number of elements along each axis");

      double gamma = 0.1;
      Sundance::setOption("gamma", gamma, "perturbation to target");

      int quadOrder = 2;
      Sundance::setOption("quadOrder", quadOrder, "Order for the Gaussian Quadrature rule");

      Sundance::init(&argc, &argv);

      VectorType<double> epetraVecType = new EpetraVectorType();

      // Define t and its derivative
      Expr dt = new Derivative(0);
      Expr t = new CoordExpr(0);

      // Create the matrix A = {{-2,1},{2,-2}}
      Expr A = List(List(-1.,0.),List(-1.,-2.));
      Expr At = List(List(-1.,-1.),List(0.,-2.));

      // Define b(t) from the MMS
      Expr b = List((2*t)/Power(E,t), Power(E,-2*t) + (0.5 + Power(t,2))/Power(E,t));
      Expr xTarget = List((0.5 + Power(t,2))/Power(E,t) + (-0.5 + t)*t*gamma,
			  (0.3333333333333333 + t)/Power(E,2*t) + 
			  ((-0.75 + t)*gamma)/5.);
      Array<double> alphaExact = Teuchos::tuple(
						0.49999999999999994 + 0.05503059750329588*gamma,
						0.3333333333333333 - 0.1272325603758784*gamma);

      // Define the mesh
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, nx, meshType);
      Mesh mesh = mesher.getMesh();

      // Set up appropriate cell filters
      CellFilter interior = new MaximalCellFilter();
      CellFilter verts = new DimensionalCellFilter(0);
      CellFilter left = verts.coordSubset(0,0.0);
      CellFilter right = verts.coordSubset(0,1.0);

      // The basis for the discrete space is tied to the dimension of the mesh
      // Even though x \n R^2, t \in R, and thus bas is 1D
      BasisFamily bas = new Lagrange(1);
      Expr x = List(new UnknownFunction(bas, "x1"), new UnknownFunction(bas, "x2"));
      Expr lambda = List(new UnknownFunction(bas), new UnknownFunction(bas));
      Expr alpha = List(new UnknownFunction(bas), new UnknownFunction(bas));

      Expr xHat = List(new TestFunction(bas), new TestFunction(bas));
      Expr lambdaHat = List(new TestFunction(bas), new TestFunction(bas));
      Expr alphaHat = List(new TestFunction(bas), new TestFunction(bas));

      QuadratureFamily quad = new GaussianQuadrature(quadOrder);


      /* state equation & BCs */
      Expr stateEqn = Integral(interior, lambdaHat*(dt*x - (A*x) - b), quad);
      Expr stateBC = EssentialBC(left, lambdaHat*(x-alpha), quad);

      /* adjoint equation & BCs, derived by hand */
      Expr adjointEqn = Integral(interior, xHat*(x-xTarget) -  xHat*(dt*lambda+(At*lambda)), quad);
      Expr adjointBC = EssentialBC(right, xHat*lambda, quad);


      /* design equation and BC */
      /* -- the (alpha')*(alphaHat') term enforces constancy of alpha in time */
      Expr designEqn = Integral(interior, (dt*alpha)*(dt*alphaHat), quad);
      Expr designBC = EssentialBC(left, alphaHat*lambda, quad);

      /* Combine the equations and BCs to form the KKT system  */
      Expr eqn = stateEqn + adjointEqn + designEqn;
      Expr bc = stateBC + adjointBC + designBC;

      /* create the LP */
      LinearProblem LP(mesh, eqn, bc, List(lambdaHat, xHat, alphaHat),
		       List(x, lambda, alpha), epetraVecType);

      ParameterList params;
      LinearSolver<double> solver = new AmesosSolver(params);
	  
      Expr soln = LP.solve(solver);

      FieldWriter writer = new DSVWriter("2DLinearOpt.dat");
      writer.addMesh(mesh);
      writer.addField("x[0]", new ExprFieldWrapper(soln[0]));
      writer.addField("x[1]", new ExprFieldWrapper(soln[1]));
      writer.addField("lambda[0]", new ExprFieldWrapper(soln[2]));
      writer.addField("lambda[1]", new ExprFieldWrapper(soln[3]));
      writer.addField("alpha[0]", new ExprFieldWrapper(soln[4]));
      writer.addField("alpha[1]", new ExprFieldWrapper(soln[5]));
      writer.write();
   
      
      Array<double> alphaNum = Teuchos::tuple(
					      L2Norm(mesh, left, soln[0], quad),
					      L2Norm(mesh, left, soln[1], quad)
					      );
      for (int j=0; j<2; j++)
	{
	  Tabs tab1;
	  Out::os() << tab1 << "Alpha[" << j << "]: exact="
		    << alphaExact[j]
		    << ", numerical " << alphaNum[j]
		    << ", error=" << fabs(alphaExact[j] - alphaNum[j])
		    << endl;
	}
    }
  catch(std::exception& ex)
    {
      cerr << "main() caught exception: " << ex.what() << endl;
    }
  Sundance::finalize();
}

    

    
    
