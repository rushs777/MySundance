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

      int quadOrder = 2;
      Sundance::setOption("quadOrder", quadOrder, "Order for the Gaussian Quadrature rule");

      Sundance::init(&argc, &argv);

      VectorType<double> epetraVecType = new EpetraVectorType();

      // Define t and its derivative
      Expr dt = new Derivative(0);
      Expr t = new CoordExpr(0);

      // Create the matrix A = {{1,2},{3,4}}
      Expr A = List(List(1.0, 2.0), List(3.0, 4.0));
      Expr At = List(List(1.0, 3.0), List(2.0, 4.0));

      // Define xTarget and  b(t) from the MMS
      //Expr b = List( -0.9747447707505879/Power(E,t),
      //		     -2.4368619268764697/Power(E,t) - 0.9391058557936693*Power(E,t));
      //      Expr xTarget = List(t, 1.0/3.0);
      Expr xTarget = List(0.5+t, 3*(1.0-t));
      Expr b = List(0. - 5.105508766453743/Power(E,t),
   -12.763771916134356/Power(E,t) - 
		    1.7459301114595882*Power(E,t));

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
      Expr x = List(new UnknownFunction(bas), new UnknownFunction(bas));
      Expr lambda = List(new UnknownFunction(bas), new UnknownFunction(bas));
      Expr alpha = List(new UnknownFunction(bas), new UnknownFunction(bas));

      Expr xHat = List(new TestFunction(bas), new TestFunction(bas));
      Expr lambdaHat = List(new TestFunction(bas), new TestFunction(bas));
      Expr alphaHat = List(new TestFunction(bas), new TestFunction(bas));

      QuadratureFamily quad = new GaussianQuadrature(quadOrder);


      /* state equation & BCs */
      Expr stateEqn = Integral(interior, lambdaHat*(dt*x - A*x - b), quad);
      //Expr stateEqn = Integral(interior, lambdaHat*(-dt*x + A*x + b), quad);
      cout << "Did stateEqn" << endl;
      Expr stateBC = EssentialBC(left, lambdaHat*(x-alpha), quad);
      cout << "Did stateEqn BC" << endl;

      /* adjoint equation & BCs, derived by hand */
      Expr adjointEqn = Integral(interior, xHat*(x-xTarget) - xHat*(dt*lambda+At*lambda), quad);
      //Expr adjointEqn = Integral(interior, xHat*(x-xTarget) + xHat*(dt*lambda+At*lambda), quad);
      cout << "Did adjointEqn" << endl;
      Expr adjointBC = EssentialBC(right, xHat*lambda, quad);
      cout << "Did adjointBC" << endl;
      
      /* Here's the hack: the design variable is only defined at an initial point. We
       * extend it to the entire time interval. To give an equation for alpha, let 
       * n be a positive integer and add
       * Integrate[eps/2 pow(t,n)*alpha*alpha, {t,0,1}] 
       * to the objective function. The factor
       * of pow(t,n)  is so that the value of alpha(0) doesn't change the value of 
       * the integral for alpha in L^2.. However, since alpha is discretized with PWLL,
       * there is a small O(h^2) contribution to the integral from this extra integral.
       * Since O(h^2) is the same order as the discretization error in the solution, and 
       * since we can multiply the "hack" term by an arbitrarily small epsilon, the
       * hack doesn't harm accuracy. 
       */
      /* Design eqn and BCs, derived by hand */
      double eps = 1.0e-4;
      Expr designEqn = Integral(interior, alphaHat*(eps*t*t*alpha), quad);
      cout << "Did designEqn" << endl;
      Expr designBC = EssentialBC(left, -alphaHat*lambda, quad);
      //Expr designBC = EssentialBC(left, alphaHat*lambda, quad);
      cout << "Did designBC" << endl;

      Expr eqn = stateEqn + adjointEqn + designEqn;
      Expr bc = stateBC + adjointBC + designBC;

      cout << "Defining the linear problem" << endl;
      LinearProblem LP(mesh, eqn, bc, List(lambdaHat, xHat, alphaHat),
		       List(x, lambda, alpha), epetraVecType);
      cout << "Linear problem defined" << endl;

      ParameterList params;
      LinearSolver<double> solver = new AmesosSolver(params);

      cout << "Right before solve" << endl;
      Expr soln = LP.solve(solver);
      cout << "Right after solve" << endl;

      FieldWriter writer = new DSVWriter("2DLinearOpt.dat");
      writer.addMesh(mesh);
      writer.addField("x[0]", new ExprFieldWrapper(soln[0]));
      writer.addField("x[1]", new ExprFieldWrapper(soln[1]));
      writer.addField("lambda[0]", new ExprFieldWrapper(soln[2]));
      writer.addField("lambda[1]", new ExprFieldWrapper(soln[3]));
      writer.addField("alpha[0]", new ExprFieldWrapper(soln[4]));
      writer.addField("alpha[1]", new ExprFieldWrapper(soln[5]));
      writer.write();




      
      
    }
  catch(std::exception& ex)
    {
      cerr << "main() caught exception: " << ex.what() << endl;
    }
  Sundance::finalize();
}

    

    
    
