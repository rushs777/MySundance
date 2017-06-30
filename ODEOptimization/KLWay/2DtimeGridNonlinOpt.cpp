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
      
      int nx = 512;
      Sundance::setOption("nx", nx, "Number of elements along each axis");

      int quadOrder = 4;
      Sundance::setOption("quadOrder", quadOrder, "Order for the Gaussian Quadrature rule");

      Sundance::init(&argc, &argv);

      VectorType<double> vecType = new EpetraVectorType();

      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, nx, meshType);
      Mesh mesh = mesher.getMesh();

      CellFilter interior = new MaximalCellFilter();
      CellFilter verts = new DimensionalCellFilter(0);
      CellFilter left = verts.coordSubset(0,0.0);
      CellFilter right = verts.coordSubset(0,1.0);

      // Define t and its derivative
      Expr dt = new Derivative(0);
      Expr t = new CoordExpr(0);

      // Create the matrix A = {{1,2},{3,4}}
      Expr A = List(List(1.0, 2.0), List(3.0, 4.0));
      Expr At = List(List(1.0, 3.0), List(2.0, 4.0));

      // Create the tensor T = { {{1,4},{0,1}}, {{3,6},{1,4}} }
      // To do T(x,x): T*x*x
      Expr T = List( List(List(1.0,4.0), List(0.0,1.0)), List(List(3.0,6.0), List(1.0,4.0)) );
      Expr Tt = List( List(List(1.0,0.0), List(4.0,1.0)), List(List(3.0,1.0), List(6.0,4.0)) );


      /* This is the new way to check operations. Create an UnknownFunction(bas, "alias")
      Expr x1 = new UnknownFunction(bas,"x1");
      Expr x2 = new UnknownFunction(bas,"x2");
      Expr test = List(x1, x2);
      cout << "Here is T*test: " << T*test << endl;
      cout << "Here is test*T*test: " << T*test*test << endl;
      */

      cout << "Here is T[0]: " << T[0] << endl;
      cout << "Here is Tt[0]: " << Tt[0] << endl;
      cout << "Here is T[0][0]: " << T[0][0] << endl;
      cout << "Here is T[0][1]: " << T[0][1] << endl;

      // Define xTarget and  b(t) from the MMS
      Expr b = List( - 0.9747447707505879/Power(E,t) - 0.31303528526455643*Power(E,t)*(0. + 0.31303528526455643*Power(E,t)) - (0.4873723853752939* (0.4873723853752939/Power(E,t) + 1.2521411410582257*Power(E,t)))/Power(E,t), -2.4368619268764697/Power(E,t) - 0.9391058557936693*Power(E,t) -  0.31303528526455643* Power(E,t)* (0.4873723853752939/Power(E,t) + 0.9391058557936693*Power(E,t)) - (0.4873723853752939*(1.9494895415011757/Power(E,t) + 1.8782117115873387*Power(E,t)))/Power(E,t));
      Expr xTarget = List(t, 1.0/3.0);

      // The basis for the discrete space is tied to the dimension of the mesh
      // Even though x \n R^2, t \in R, and thus bas is 1D
      BasisFamily bas = new Lagrange(1);
      Expr x = List(new UnknownFunction(bas), new UnknownFunction(bas));
      Expr lambda = List(new UnknownFunction(bas), new UnknownFunction(bas));
      Expr alpha = List(new UnknownFunction(bas), new UnknownFunction(bas));

      Expr xHat = List(new TestFunction(bas), new TestFunction(bas));
      Expr lambdaHat = List(new TestFunction(bas), new TestFunction(bas));
      Expr alphaHat = List(new TestFunction(bas), new TestFunction(bas));

      QuadratureFamily quad = new GaussianQuadrature(quadOrder);

      /* state eqn & bcs, derived by hand */
      Expr stateEqn = Integral(interior, lambdaHat*(dt*x - T*x*x - A*x -b), quad);
      Expr stateBC = EssentialBC(left, lambdaHat*(x-alpha), quad);

      /* adjoint eqn & bcs, derived by hand */
      Expr DerivT = List( List( (T[0]*x + Tt[0]*x)[0], (T[0]*x + Tt[0]*x)[1] ),
			  List( (T[1]*x + Tt[1]*x)[0], (T[1]*x + Tt[1]*x)[1] ));
      Expr adjointEqn = Integral(interior, xHat*(x-xTarget) - xHat*(dt*lambda+At*lambda - DerivT*lambda), quad);
      Expr adjointBC = EssentialBC(right, xHat*lambda, quad);

      /* design eqn & bcs */
      /* See the note in timeGridOpt.cpp on the hack used in the design equation*/
      double eps = 1.0e-4;
      Expr designEqn = Integral(interior, alphaHat*(eps*t*t*alpha), quad);
      Expr designBC = EssentialBC(left, -alphaHat*lambda, quad);

      Expr eqn = stateEqn + adjointEqn + designEqn;
      Expr bc = stateBC + adjointBC + designBC;

      cout << "Finished setting up equations " << endl;
      
      DiscreteSpace ds(mesh, List(bas, bas, bas), vecType);
      Expr U0 = new DiscreteFunction(ds, 0.0);
      
      NonlinearProblem NLP(mesh, eqn, bc, List(lambdaHat, xHat, alphaHat),
			   List(x, lambda, alpha), U0, vecType);
      
      NonlinearSolver<double> solver 
	= NonlinearSolverBuilder::createSolver("playa-newton-amesos.xml");

      SolverState<double> state = NLP.solve(solver);
    
      TEUCHOS_TEST_FOR_EXCEPTION(state.finalState() != SolveConverged,
				 std::runtime_error,
				 "Nonlinear solve failed to converge: message="
				 << state.finalMsg());

      //Expr f = Integral(interior, 0.5*(U0[0]-t)*(U0[0]-t),quad);
      //Vector<double> fVec = getDiscreteFunctionVector(f);
      //Expr f_discrete = new DiscreteFunction(ds, fVec);      
      
      FieldWriter writer = new DSVWriter("2DNonlin-opt.dat");
      writer.addMesh(mesh);
      writer.addField("x", new ExprFieldWrapper(U0[0]));
      writer.addField("lambda", new ExprFieldWrapper(U0[1]));
      writer.addField("alpha", new ExprFieldWrapper(U0[2]));
      //writer.addField("f", new ExprFieldWrapper(f_discrete));
      writer.write();
      

      
    }
  catch(std::exception& ex)
    {
      cerr << "main() caught exception: " << ex.what() << endl;
    }
  Sundance::finalize();
}

    

    
    
