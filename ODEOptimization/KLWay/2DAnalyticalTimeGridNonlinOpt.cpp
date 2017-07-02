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

      Expr xTarget = List(1/(2.*Power(E,t)) + ((-0.5 + t)*gamma)/4.,
			  Power(E,-0.5 + 1/(2.*Power(E,t)) + t)/3. + 
			  (1 + t)*gamma);
      Expr b = List(0.0,0.0);

      Expr A = List( List(-1.0,0.0), List(0.0,1.0) );
      Expr At = A; //A is symmetric
      Expr T = List( List( List(0.0,0.0), List(0.0,0.0) ), List( List(0.0,-1.0), List(0.0,0.0) ) );
      
      /* "exact solution" produced by Mathematica's FindMinimum. Accurate to
       * about 1.0e-6. */
      Array<double> alphaExact = Teuchos::tuple(
						0.5000000208134503 - 0.027365139770605805*gamma + 
						0.00960551132676607*gamma*gamma + 
						0.0025004038373395636*gamma*gamma*gamma,
						0.3333333219069475 + 1.0656329517752654*gamma - 
						0.01145702377095341*gamma*gamma + 
						0.005347348624291412*gamma*gamma*gamma
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
      /*Expr stateEqn = Integral(interior,
			       lambdaHat[0]*(dt*x[0] + x[0])
			       + lambdaHat[1]*(dt*x[1] - x[1] + x[0]*x[1])
			       , quad);*/
      Expr stateEqn = Integral(interior, lambdaHat*(dt*x - A*x - T*x*x -b), quad);
      Expr stateBC = EssentialBC(left, lambdaHat*(x-alpha), quad);

      /* adjoint equation & BCs, derived by hand */
      /*Expr adjointEqn = Integral(interior, xHat*(x-xTarget) 
				 - xHat[0]*(dt*lambda[0] - lambda[0] - lambda[1]*x[1])
				 - xHat[1]*(dt*lambda[1] + lambda[1] - lambda[1]*x[0]),
				 quad);*/
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

      
      FieldWriter writer = new DSVWriter("2DNonlinearOpt.dat");
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

    

    
    
