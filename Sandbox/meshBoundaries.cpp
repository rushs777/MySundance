#include "Sundance.hpp"
#include "meshAddOns.hpp"


using std::cout;
using std::endl;



/**
 * This class is designed to give a CELL_PREDICATE test that
 * returns true for all Points that do not have an x-coordinate 
 * equal to dim0_max_
 */
class outflowEdgeTest: public CellPredicateFunctorBase,
		       public Playa::Handleable<CellPredicateFunctorBase>
{
public:
  outflowEdgeTest(const double& dim0_max) : CellPredicateFunctorBase("outflowEdgeTest"), dim0_max_(dim0_max) {}
  virtual ~outflowEdgeTest() {}
  virtual bool operator()(const Point& x) const {return fabs(x[0] - dim0_max_) > 1.0e-10;}
  GET_RCP(CellPredicateFunctorBase);

private:
  double dim0_max_;
};




int main(int argc, char *argv[])
{
  int nx = 3;
  Sundance::setOption("nx",nx,"Number of cells along each axis");

  int meshVerbosity = 0;
  Sundance::setOption("meshVerbosity",meshVerbosity,"Level of verbosity during mesh creation");

  Sundance::init(&argc, &argv);

  
  // Define the boundaries of the mesh
  double xmin = 0.0;
  double xmax = 1.0;
  double ymin = 0.0;
  double ymax = 1.0;

  // Create the 2D mesh
  MeshType meshType = new Sundance::BasicSimplicialMeshType();
  MeshSource mesher = new Sundance::PartitionedRectangleMesher(xmin, xmax, nx, 1,
							       ymin, ymax, nx, 1,
							       meshType, meshVerbosity);
  Mesh mesh = mesher.getMesh();

  displayVertices(mesh);

  Array<double> boundaries = getBoundaries(mesh);
  cout << "The boundaries of the mesh are: ["
       << boundaries[0] << "," << boundaries[1] << "] x ["
       << boundaries[2] << "," << boundaries[3] << "]" << endl;

  CellFilter bdryFilter = new BoundaryCellFilter();
  CellFilter GammaStar = bdryFilter.subset( new outflowEdgeTest(boundaries[1]) );

  
  cout << mesher.description() << endl;

  return 0;
}
