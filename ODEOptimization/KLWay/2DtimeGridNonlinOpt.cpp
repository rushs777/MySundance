#include "Sundance.hpp"
#include "PlayaAmesosSolver.hpp"
#include "PlayaDenseSerialMatrix.hpp"
#include "PlayaSerialVectorType.hpp"

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
      Expr T = List(List(List(0.05,0.1),List(0.15,0.2)),
		    List(List(0.05,0.05),List(0.,0.15)));

      // Define b(t) from the MMS
      Expr b = List((40*Power(E,3*t)*t - 
		     4*Power(0.3333333333333333 + t,2) - 
		     5*Power(E,t)*(0.3333333333333333 + t)*
		     (0.5 + Power(t,2)) - 
		     Power(E,2*t)*Power(0.5 + Power(t,2),2))/
		    (20.*Power(E,4*t)),
		    (-4*Power(1 + 3*t,2) + 
		     120*Power(E,3*t)*(1 + 2*Power(t,2)) - 
		     2*Power(E,t)*(1 + 3*t)*(1 + 2*Power(t,2)) - 
		     3*Power(E,2*t)*(-79 + 
				     4*(Power(t,2) + Power(t,4))))/
		    (240.*Power(E,4*t)));
      Expr xTarget = List((0.5 + Power(t,2))/Power(E,t) + (-0.5 + t)*t*gamma,
			  (0.3333333333333333 + t)/Power(E,2*t) + 
			  ((-0.75 + t)*gamma)/5.);
      /* "exact solution" produce by Mathematica's FindMinimum. Accurate to
       * about 1.0e-6. */
      Array<double> alphaExact = Teuchos::tuple(
						0.5000000289639744 + gamma*
						(0.064717637614702 + 
						 (-0.00044913462946438883 - 
						  0.00001613396940867341*gamma)*gamma),
						0.3333333245793077 + gamma*
						(-0.10739738478174649 + 
						 (-0.0015852504354095308 - 
						  0.000019917356023471516*gamma)*gamma)
						);

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

      Array<BasisFamily> basArray(6);
      for (int i=0; i<basArray.size(); i++) basArray[i]=bas;
      DiscreteSpace ds(mesh, basArray, epetraVecType);
      Expr U0 = new DiscreteFunction(ds, 0.0);

      QuadratureFamily quad = new GaussianQuadrature(quadOrder);


      /* state equation & BCs */
      /*Expr stateEqn = Integral(interior, lambdaHat*(dt*x - (A*x) - b), quad);
      for (int i=0; i<2; i++)
	{
	  for (int j=0; j<2; j++)
	    {
	      for (int k=0; k<2; k++)
		{
		  stateEqn = stateEqn + Integral(interior,
						-lambdaHat[i]*T[i][j][k]*x[j]*x[k], quad);
		}
	    }
	    }*/
      Expr stateEqn = Integral(interior, lambdaHat*(dt*x - A*x - T*x*x -b), quad);
      Expr stateBC = EssentialBC(left, lambdaHat*(x-alpha), quad);

      /* adjoint equation & BCs, derived by hand */
      /*Expr adjointEqn = Integral(interior, xHat*(x-xTarget) -  xHat*(dt*lambda+(At*lambda)), quad);
      for (int i=0; i<2; i++)
	{
	  for (int j=0; j<2; j++)
	    {
	      for (int k=0; k<2; k++)
		{
		  adjointEqn = adjointEqn + Integral(interior,
						     lambda[i]*(T[i][j][k]+T[i][k][j])*x[j]*xHat[k], quad);
		}
	    }
	    }*/
      Expr adjointEqn = Integral(interior, xHat*(x-xTarget) - xHat*(dt*lambda+At*lambda) - lambda*(T*xHat*x + T*x*xHat), quad);
      Expr adjointBC = EssentialBC(right, xHat*lambda, quad);

      /* design equation and BC */
      /* -- the (alpha')*(alphaHat') term enforces constancy of alpha in time */
      Expr designEqn = Integral(interior, (dt*alpha)*(dt*alphaHat), quad);
      Expr designBC = EssentialBC(left, alphaHat*lambda, quad);

      /* Combine the equations and BCs to form the KKT system  */
      Expr eqn = stateEqn + adjointEqn + designEqn;
      Expr bc = stateBC + adjointBC + designBC;

      /* create the NLP */
      NonlinearProblem NLP(mesh, eqn, bc, List(lambdaHat, xHat, alphaHat),
			   List(x, lambda, alpha), U0, epetraVecType);

      NonlinearSolver<double> solver 
        = NonlinearSolverBuilder::createSolver("playa-newton-amesos.xml");

      SolverState<double> state = NLP.solve(solver);
      

      TEUCHOS_TEST_FOR_EXCEPTION(state.finalState() != SolveConverged,
                                 std::runtime_error,
                                 "Nonlinear solve failed to converge: message="
                                 << state.finalMsg());

      
      FieldWriter writer = new DSVWriter("2DLinearOpt-.dat");
      writer.addMesh(mesh);
      writer.addField("x[0]", new ExprFieldWrapper(U0[0]));
      writer.addField("x[1]", new ExprFieldWrapper(U0[1]));
      writer.addField("lambda[0]", new ExprFieldWrapper(U0[2]));
      writer.addField("lambda[1]", new ExprFieldWrapper(U0[3]));
      writer.addField("alpha[0]", new ExprFieldWrapper(U0[4]));
      writer.addField("alpha[1]", new ExprFieldWrapper(U0[5]));
      writer.write();
      
      Array<double> alphaNum = Teuchos::tuple(
					      L2Norm(mesh, left, U0[0], quad),
					      L2Norm(mesh, left, U0[1], quad)
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

    

    
    
