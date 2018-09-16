#include "Sundance.hpp"

Array<Expr> GramSchmidt(Array<Expr> v, Mesh mesh, int quadOrder)
{
  if(v.length()<1)
    {
      Out::root() << "Array of length 0 passed to GramSchmidt" << endl;
      return v;
    }

  CellFilter interior = new MaximalCellFilter();
  QuadratureFamily quad = new GaussianQuadrature(quadOrder);

  v[0] = v[0]*(1.0/L2Norm(mesh, interior, v[0], quad));

  for(int i=1; i<v.length(); i++)
    {
      // make v[i] orthogonal
      for(int j=0; j<i; j++)
	{
	  FunctionalEvaluator IP = FunctionalEvaluator(mesh, Integral(interior, v[i]*v[j], quad));
	  v[i] = v[i] - IP.evaluate()*v[j];
	}

      // make v[i] orthonormal
      v[i] = v[i]*(1.0/L2Norm(mesh, interior, v[i], quad));
    }

  return v;
}

int main(int argc, char *argv[])
{
  try
    {
      int nx = 2;
      Sundance::setOption("nx",nx,"Number of elemens along each axis");

      Sundance::init(&argc, &argv);
      
      // Define the coordinates
      Expr x = new CoordExpr(0,"x");
      Expr y = new CoordExpr(1,"y");
      
      // Define the mesh
      MeshType meshType = new BasicSimplicialMeshType();
      double xmin = -1.0;
      double xmax = 1.0;
      double ymin = 0.0;
      double ymax = 1.0;
      MeshSource mesher = new PartitionedRectangleMesher(xmin,xmax,nx,1,ymin,ymax,nx,1,meshType,2);
      Mesh mesh = mesher.getMesh();

      MeshSource mesher1D = new PartitionedLineMesher(xmin,xmax,nx,meshType,0);
      Mesh mesh1D = mesher1D.getMesh();

      Expr f = List(x*x+y, y);

      cout << "f: " << f << endl;

      BasisFamily basisType = new Sundance::Lagrange(2);
      Array<BasisFamily> basis = List(basisType, basisType);
      Playa::VectorType<double> epetraVecType = new Playa::EpetraVectorType();
      DiscreteSpace ds(mesh, basis, epetraVecType);

      L2Projector fproj(ds, f);
      Expr projected_f = fproj.project();

      //cout << "Here is projected_f: " << projected_f << endl;

      Vector<double> fVec = getDiscreteFunctionVector(projected_f);

      //cout << "Here is fVec: " << endl << fVec << endl;
      //cout << "Number of edges plus vertices: " << mesh.numCells(0) + mesh.numCells(1) << endl;

      Expr one = 1.0;
      Array<Expr> GS_test = Teuchos::tuple(sin(x), x, cos(x));
      cout << "v = " << GS_test << endl;

      Array<Expr> GS_done = GramSchmidt(GS_test, mesh1D, 4);

      cout << "u = " << GS_done << endl;
      CellFilter interior = new MaximalCellFilter;
      QuadratureFamily quad = new GaussianQuadrature(4);

      for(int i=0; i<GS_done.length(); i++)
	{
	  for(int j=i+1; j<GS_done.length(); j++)
	    {
	      FunctionalEvaluator IP = FunctionalEvaluator(mesh, Integral(interior, GS_done[i]*GS_done[j], quad));
	      cout << "<u[" <<i << "],u[" << j << "]> = " << IP.evaluate() << endl;
	    }
	}

      /*      Array<int> array1 = tuple(1,2,3);
      Array<int> array2 = array1;
      cout << "array1 = " << array1 << endl;
      cout << "array2 = " << array2 << endl;
      array1[2] = 5;
      array2[2] = 10;
      cout << "array1 = " << array1 << endl;
      cout << "array2 = " << array2 << endl;
      */

      cout << "mesh's ordering convention: " << mesh.meshOrder() << endl;
      cout << "Number of vertices: " << mesh.numCells(0) << endl;
      cout << "number of edges: " << mesh.numCells(1) << endl;

      cout << "mesher's description: " << mesher.description() << endl;

      Expr integrand = f[0];
      cout << "integrand: " << integrand << endl;
      Expr eqn = Integral(interior,integrand,quad);
      FunctionalEvaluator IP = FunctionalEvaluator(mesh, eqn);
      cout << "Value of integral over [-1,1]x[0,1]: " << IP.evaluate() << endl;

      CellFilter slice = interior.subset(new CoordinateValueCellPredicate(0,0.5));
      Expr eqn2 = Integral(integralDomain, integrand,quad);
      FunctionalEvaluator integrator = FunctionalEvaluator(mesh, eqn2);
      cout << "Value of integral over {0.5}x[0,1]: " << integrator.evaluate() << endl;
    }
  catch(std::exception& e)
    {
      Out::root() << "Caught exception: " << e.what() << std::endl;
      return -1;
    }  
}
