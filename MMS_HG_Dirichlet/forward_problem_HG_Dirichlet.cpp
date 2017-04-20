#include "VientoPCDStepper.hpp"
#include "VientoErrorChecker.hpp"
#include "Sundance.hpp"
#include "VientoDefaultOutputManager.hpp"



using namespace Sundance;
using namespace Viento;
using Sundance::List;

// Things necessary for mathematica generated fs to be understood
#include "MathematicaConverter.hpp"

/** 
 * MMS flow with a projection method
 */

CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-2.0*M_PI) < 1.0e-10;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-2.0*M_PI) < 1.0e-10;})

CELL_PREDICATE(CornerPointTest, {return fabs(x[1]) < 1.0e-10 && fabs(x[0])<1.0e-10;})

Expr uExFunc(const Expr& t)
{
  Expr x = new CoordExpr(0);
  Expr y = new CoordExpr(1);

  return List((2*Cos(2*t)*Cos(y)*Power(Sin(x),2)*Sin(y))/3. + 
    		      (Cos(3*t)*Cos(y)*Power(Sin(2*x),2)*Sin(y))/2.,
		      (-2*Cos(2*t)*Cos(x)*Sin(x)*Power(Sin(y),2))/3. - 
		      Cos(3*t)*Cos(2*x)*Sin(2*x)*Power(Sin(y),2));
}

Expr pExFunc(const Expr& t)
{
  Expr x = new CoordExpr(0);
  Expr y = new CoordExpr(1);

  return Cos(x)*Cos(y)*Sin(t);
}

Expr qExFunc(const Expr& t)
{
  Expr x = new CoordExpr(0);
  Expr y = new CoordExpr(1);
  double nu = 1.0;

  return List((-36*Cos(y)*Sin(t)*Sin(x) + (((44 + 30*Cos(t) + 8*Cos(4*t) + 30*Cos(5*t))*
             Cos(x) + 9*((3 + 2*Cos(t) + 2*Cos(5*t))*Cos(3*x) + Cos(5*x)))*
          Power(Sin(x),3) + 9*Cos(6*t)*Cos(2*x)*Power(Sin(2*x),3))*Power(Sin(y),2) - 
      3*(8*nu*Cos(2*t)*(-1 + 2*Cos(2*x)) + 6*nu*Cos(3*t)*(-1 + 5*Cos(4*x)) + 
         8*Sin(2*t)*Power(Sin(x),2) + 9*Sin(3*t)*Power(Sin(2*x),2))*Sin(2*y))/36.,
   (6*nu*(2*Cos(2*t)*(-1 + 2*Cos(2*y))*Sin(2*x) + 
         3*Cos(3*t)*(-4 + 5*Cos(2*y))*Sin(4*x)) - 18*Cos(x)*Sin(t)*Sin(y) + 
      3*(4*Sin(2*t)*Sin(2*x) + 9*Sin(3*t)*Sin(4*x))*Power(Sin(y),2) + 
      2*Cos(y)*(Cos(2*t)*(4*Cos(2*t) - 3*Cos(3*t)*(-3 - 6*Cos(2*x) + Cos(4*x)))*
          Power(Sin(x),2) + 9*Power(Cos(3*t),2)*Power(Sin(2*x),2))*Power(Sin(y),3))/
    18.);
}





int main(int argc, char** argv)
{
  try
    {
      Time timer("total");
      timer.start();
      int nx = 32;
      Sundance::setOption("nx", nx, "grid size");

      int nSteps = 32;
      Sundance::setOption("nSteps", nSteps, "number of timesteps");

      double tFinal = 1.0;
      Sundance::setOption("tFinal", tFinal, "final time");

      int order=2;
      Sundance::setOption("order", order, "stepper order");

      int verb = 1;
      Sundance::setOption("verbosity", verb, "verbosity");
      Sundance::init(&argc, &argv);

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedRectangleMesher. */
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedRectangleMesher(0.0, 2.0*M_PI, nx, 1,
                                                         0.0, 2.0*M_PI, nx, 1,
                                                         meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter bdry = new BoundaryCellFilter();
      CellFilter top = bdry.subset( new TopPointTest() );
      CellFilter right = bdry.subset( new RightPointTest() );
      CellFilter bottom = bdry.subset( new BottomPointTest() );
      CellFilter left = bdry.subset( new LeftPointTest() );
      CellFilter nodes = new DimensionalCellFilter(0);
      CellFilter corner = nodes.subset(new CornerPointTest());

      /* Bases for velocity and pressure */
      BasisFamily vBas  = new Lagrange(2);
      BasisFamily pBas  = new Lagrange(1);

      /* Create the FlowState object */

      RCP<ProblemDescription> prob = rcp(new ProblemDescription(mesh));

      prob->addCondition(NoSlip, bdry, StateDependentFunction(uExFunc));
      prob->addCondition(BodyForce, interior, StateDependentFunction(qExFunc));
      prob->peg(corner, StateDependentFunction(pExFunc));
      prob->setNuEff(1.0);

      PCDControlParams pcdParams;
      NewtonControlParams newtonParams;
      pcdParams.FParams.method = "CG";
      PCDStepper stepper(prob, pcdParams, newtonParams, order, SemiImplicit);
      RCP<FlowState> state = stepper.state();

      double dt = tFinal/nSteps;
      double tCur = 0.0;
      
      state->initialize(tCur, dt, uExFunc, pExFunc);

      string outDir = "Results/ForwardProblem/";
   //   system( ("mkdir -p " + outDir).c_str() ); Done in DefaultOutputManager constructor
      string filename = outDir + "HG_Dirichlet_nx" + Teuchos::toString(nx)
	+ "_nt" + Teuchos::toString(nSteps);

      ErrorChecker check(state, uExFunc, pExFunc);

//Added by me to write out the snapshots
      RCP<DefaultOutputManager> output
	= rcp(new DefaultOutputManager(filename, "st", new ExodusWriterFactory()));
	
      for (int i=0; i<=nSteps; i++)
	{
	  Tabs tab1;
	  SUNDANCE_MSG1(verb, tab1<< "step " << i << " time=" << tCur);
	  if (i>0)
	    {
	      stepper.step(dt, verb);
	    }
	  Tabs tab2;
	  SUNDANCE_MSG1(verb, tab2<< "CFL=" << state->minHOverV() );
	  FieldWriter w = new VTKWriter(filename + "/" + Teuchos::toString(i));
	  check.write(w);
	  output->write(i, stepper.state(), 1);
	}

      timer.stop();
      Out::root() << "nx=" << nx << " nt=" << nSteps << " uErr=" << check.uErrL2()
		  << " pErr=" << check.pErrL2() << " runtime=" << timer.totalElapsedTime() << endl;

    }
  catch(std::exception& e)
    {
      Sundance::handleException(e);
    }
  Sundance::finalize(); return Sundance::testStatus(); 
}
