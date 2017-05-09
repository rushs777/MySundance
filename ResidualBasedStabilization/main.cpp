#include "Sundance.hpp"

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
      double xmin = 0.0;
      double xmax = 1.0;
      double ymin = 0.0;
      double ymax = 1.0;
      MeshSource mesher = new PartitionedRectangleMesher(xmin,xmax,nx,1,ymin,ymax,nx,1,meshType,0);
      Mesh mesh = mesher.getMesh();

      Expr f = List(x*x, y);

      cout << "Hi, here is f: " << f << endl;

      BasisFamily basisType = new Sundance::Lagrange(2);
      Array<BasisFamily> basis = List(basisType, basisType);
      Playa::VectorType<double> epetraVecType = new Playa::EpetraVectorType();
      DiscreteSpace ds(mesh, basis, epetraVecType);

      L2Projector fproj(ds, f);
      Expr projected_f = fproj.project();

      cout << "Here is projected_f: " << projected_f << endl;

      Vector<double> fVec = getDiscreteFunctionVector(projected_f);

      cout << "Here is fVec: " << endl << fVec << endl;

      cout << "Number of edges plus vertices: " << mesh.numCells(0) + mesh.numCells(1) << endl;
    }
  catch(std::exception& e)
    {
      Out::root() << "Caught exception: " << e.what() << std::endl;
      return -1;
    }  
}
