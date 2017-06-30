#include "Sundance.hpp"
#include "PlayaAmesosSolver.hpp"


int main(int argc, char *argv[]) 
{
  try
    {
      int verb = 0;      
      Sundance::setOption("verbosity", verb, "verbosity level");
      
      int nx = 512;
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

      QuadratureFamily quad = new GaussianQuadrature(4);

      /* state eqn & bcs, derived by hand */
      Expr stateEqn = Integral(interior, lambdaHat*(dt*x - x*(1.0-x)), quad);
      Expr stateBC = EssentialBC(left, lambdaHat*(x-alpha), quad);

      /* adjoint eqn & bcs, derived by hand */
      Expr adjointEqn = Integral(interior, xHat*(x-t) - xHat*(dt*lambda+(1.0-2.0*x)*lambda), quad);
      Expr adjointBC = EssentialBC(right, xHat*lambda, quad);

      /* design eqn & bcs */
      /* See the note in timeGridOpt.cpp on the hack used in the design equation*/
      double eps = 1.0e-4;
      Expr designEqn = Integral(interior, eps*alphaHat*t*t*alpha, quad);
      Expr designBC = EssentialBC(left, -alphaHat*lambda, quad);

      Expr eqn = stateEqn + adjointEqn + designEqn;
      Expr bc = stateBC + adjointBC + designBC;

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
      
      FieldWriter writer = new DSVWriter("nonlin-opt.dat");
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

    

    
    
