#include "Sundance.hpp"
#include "PlayaAmesosSolver.hpp"


int main(int argc, char *argv[]) 
{
  try
    {
      int verb = 0;      
      Sundance::setOption("verbosity", verb, "verbosity level");
      
      int nx = 128;
      Sundance::setOption("nx", nx, "Number of elements along each axis");

      Sundance::init(&argc, &argv);

      VectorType<double> vecType = new EpetraVectorType();

      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, nx, meshType);
      Mesh mesh = mesher.getMesh();

      CellFilter interior = new MaximalCellFilter();
      CellFilter verts = new DimensionalCellFilter(0);
      CellFilter left = verts.coordSubset(0,0.0);
      CellFilter right = verts.coordSubset(0,1.0);

      Expr dt = new Derivative(0);
      Expr t = new CoordExpr(0);

      BasisFamily bas = new Lagrange(1);
      Expr x = new UnknownFunction(bas);
      Expr lambda = new UnknownFunction(bas);
      Expr alpha = new UnknownFunction(bas);

      Expr xHat = new TestFunction(bas);
      Expr lambdaHat = new TestFunction(bas);
      Expr alphaHat = new TestFunction(bas);

      QuadratureFamily quad = new GaussianQuadrature(2);

      /* state equation & BCs */
      Expr stateEqn = Integral(interior, lambdaHat*(dt*x - x), quad);
      Expr stateBC = EssentialBC(left, lambdaHat*(x-alpha), quad);

      /* adjoint equation & BCs, derived by hand */
      Expr adjointEqn = Integral(interior, xHat*(x-t) - xHat*(dt*lambda+lambda), quad);
      Expr adjointBC = EssentialBC(right, xHat*lambda, quad);

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
       * have doesn't harm accuracy. 
       */
      /* Design eqn and BCs, derived by hand */
      double eps = 1.0e-4;
      Expr designEqn = Integral(interior, eps*alphaHat*t*t*alpha, quad);
      Expr designBC = EssentialBC(left, -alphaHat*lambda, quad);

      Expr eqn = stateEqn + adjointEqn + designEqn;
      Expr bc = stateBC + adjointBC + designBC;
      
      LinearProblem LP(mesh, eqn, bc, List(lambdaHat, xHat, alphaHat),
		       List(x, lambda, alpha), vecType);
      

      ParameterList params;
      LinearSolver<double> solver = new AmesosSolver(params);

      Expr soln = LP.solve(solver);

      FieldWriter writer = new DSVWriter("opt.dat");
      writer.addMesh(mesh);
      writer.addField("x", new ExprFieldWrapper(soln[0]));
      writer.addField("lambda", new ExprFieldWrapper(soln[1]));
      writer.addField("alpha", new ExprFieldWrapper(soln[2]));
      writer.write();
      

      
    }
  catch(std::exception& ex)
    {
      cerr << "main() caught exception: " << ex.what() << endl;
    }
  Sundance::finalize();
}

    

    
    
