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
 * MMS Transient Channel with a projection method
 * Domain is [0,1]x[0,1]
 */

CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-1.0) < 1.0e-10;})

CELL_PREDICATE(CornerPointTest, {return fabs(x[1]) < 1.0e-10 && fabs(x[0])<1.0e-10;})

Expr uExFunc(const Expr& t)
{
  Expr x = new CoordExpr(0);
  Expr y = new CoordExpr(1);

  return List(1 - (Pi*(-1 + Power(x,2))*
       (8*Sin((Pi*t)/2.)*
          Power(Sin(Pi*x),2)*
          Sin(4*Pi*y) + 
         45*Sin(Pi*t)*
          Power(Sin(2*Pi*x),2)*
          Sin(6*Pi*y)))/200.,
   (4*Sin((Pi*t)/2.)*Sin(Pi*x)*
       (Pi*(-1 + Power(x,2))*
          Cos(Pi*x) + x*Sin(Pi*x))*
       Power(Sin(2*Pi*y),2) + 
      15*Sin(Pi*t)*Sin(2*Pi*x)*
       (2*Pi*(-1 + Power(x,2))*
          Cos(2*Pi*x) + x*Sin(2*Pi*x)
         )*Power(Sin(3*Pi*y),2))/100.);
}

Expr pExFunc(const Expr& t)
{
  Expr x = new CoordExpr(0);
  Expr y = new CoordExpr(1);

  return 0.0;
}

Expr qExFunc(const Expr& t)
{
  Expr x = new CoordExpr(0);
  Expr y = new CoordExpr(1);
  double R = 0.0;

  return List((-(Power(Pi,2)*
         (-1 + Power(x,2))*
         (4*Cos((Pi*t)/2.)*
            Power(Sin(Pi*x),2)*
            Sin(4*Pi*y) + 
           45*Cos(Pi*t)*
            Power(Sin(2*Pi*x),2)*
            Sin(6*Pi*y))) - 
      4*Power(Pi,3)*
       (-1 + Power(x,2))*
       (32*Sin((Pi*t)/2.)*
          Power(Sin(Pi*x),2)*
          Sin(4*Pi*y) + 
         405*Sin(Pi*t)*
          Power(Sin(2*Pi*x),2)*
          Sin(6*Pi*y)) + 
      Pi*(16*Power(Pi,2)*
          (-1 + Power(x,2))*
          Cos(2*Pi*x)*Sin((Pi*t)/2.)*
          Sin(4*Pi*y) + 
         16*Sin((Pi*t)/2.)*Sin(Pi*x)*
          (4*Pi*x*Cos(Pi*x) + 
            Sin(Pi*x))*Sin(4*Pi*y) + 
         45*Sin(Pi*t)*
          (1 + 
            (-1 + 
              8*Power(Pi,2)*
              (-1 + Power(x,2)))*
             Cos(4*Pi*x) + 
            8*Pi*x*Sin(4*Pi*x))*
          Sin(6*Pi*y)) + 
      (Pi*R*(-((-1 + Power(x,2))*
              (32*Pi*Cos(4*Pi*y)*
              Sin((Pi*t)/2.)*
              Power(Sin(Pi*x),2) + 
              270*Pi*Cos(6*Pi*y)*
              Sin(Pi*t)*
              Power(Sin(2*Pi*x),2))*
              (4*Sin((Pi*t)/2.)*
              Sin(Pi*x)*
              (Pi*(-1 + Power(x,2))*
              Cos(Pi*x) + x*Sin(Pi*x)
              )*Power(Sin(2*Pi*y),2)\
              + 15*Sin(Pi*t)*
              Sin(2*Pi*x)*
              (2*Pi*
              (-1 + Power(x,2))*
              Cos(2*Pi*x) + 
              x*Sin(2*Pi*x))*
              Power(Sin(3*Pi*y),2)))\
            + (8*Sin((Pi*t)/2.)*
              Sin(Pi*x)*
              (Pi*(-1 + Power(x,2))*
              Cos(Pi*x) + x*Sin(Pi*x)
              )*Sin(4*Pi*y) + 
              45*Sin(Pi*t)*
              Sin(2*Pi*x)*
              (2*Pi*
              (-1 + Power(x,2))*
              Cos(2*Pi*x) + 
              x*Sin(2*Pi*x))*
              Sin(6*Pi*y))*
            (8*Pi*(-1 + Power(x,2))*
              Sin((Pi*t)/2.)*
              Power(Sin(Pi*x),2)*
              Sin(4*Pi*y) + 
              5*
              (-40 + 
              9*Pi*(-1 + Power(x,2))*
              Sin(Pi*t)*
              Power(Sin(2*Pi*x),2)*
              Sin(6*Pi*y)))))/100.)/
    200.,(-32*Power(Pi,2)*
       Cos(4*Pi*y)*Sin((Pi*t)/2.)*
       Sin(Pi*x)*
       (Pi*(-1 + Power(x,2))*
          Cos(Pi*x) + x*Sin(Pi*x)) - 
      270*Power(Pi,2)*Cos(6*Pi*y)*
       Sin(Pi*t)*Sin(2*Pi*x)*
       (2*Pi*(-1 + Power(x,2))*
          Cos(2*Pi*x) + x*Sin(2*Pi*x)
         ) + 2*Pi*Cos((Pi*t)/2.)*
       Sin(Pi*x)*
       (Pi*(-1 + Power(x,2))*
          Cos(Pi*x) + x*Sin(Pi*x))*
       Power(Sin(2*Pi*y),2) + 
      4*Power(Pi,2)*Sin((Pi*t)/2.)*
       Sin(Pi*x)*
       (Pi*(-1 + Power(x,2))*
          Cos(Pi*x) + x*Sin(Pi*x))*
       Power(Sin(2*Pi*y),2) - 
      4*Pi*Sin((Pi*t)/2.)*Sin(Pi*x)*
       ((4 - Power(Pi,2)*
             (-1 + Power(x,2)))*
          Cos(Pi*x) - 
         5*Pi*x*Sin(Pi*x))*
       Power(Sin(2*Pi*y),2) - 
      8*Pi*Cos(Pi*x)*Sin((Pi*t)/2.)*
       (3*Pi*x*Cos(Pi*x) + 
         (1 - Power(Pi,2)*
             (-1 + Power(x,2)))*
          Sin(Pi*x))*
       Power(Sin(2*Pi*y),2) + 
      15*Pi*Cos(Pi*t)*Sin(2*Pi*x)*
       (2*Pi*(-1 + Power(x,2))*
          Cos(2*Pi*x) + x*Sin(2*Pi*x)
         )*Power(Sin(3*Pi*y),2) + 
      60*Power(Pi,2)*Sin(Pi*t)*
       Sin(2*Pi*x)*
       (2*Pi*(-1 + Power(x,2))*
          Cos(2*Pi*x) + x*Sin(2*Pi*x)
         )*Power(Sin(3*Pi*y),2) + 
      60*Pi*Sin(Pi*t)*Sin(2*Pi*x)*
       (2*(-1 + 
            Power(Pi,2)*
             (-1 + Power(x,2)))*
          Cos(2*Pi*x) + 
         5*Pi*x*Sin(2*Pi*x))*
       Power(Sin(3*Pi*y),2) - 
      60*Pi*Cos(2*Pi*x)*Sin(Pi*t)*
       (6*Pi*x*Cos(2*Pi*x) + 
         (1 - 4*Power(Pi,2)*
             (-1 + Power(x,2)))*
          Sin(2*Pi*x))*
       Power(Sin(3*Pi*y),2) + 
      (R*((4*Sin((Pi*t)/2.)*
              Sin(Pi*x)*
              (Pi*(-1 + Power(x,2))*
              Cos(Pi*x) + x*Sin(Pi*x)
              )*Power(Sin(2*Pi*y),2)\
              + 15*Sin(Pi*t)*
              Sin(2*Pi*x)*
              (2*Pi*
              (-1 + Power(x,2))*
              Cos(2*Pi*x) + 
              x*Sin(2*Pi*x))*
              Power(Sin(3*Pi*y),2))*
            (8*Pi*Sin((Pi*t)/2.)*
              Sin(Pi*x)*
              (Pi*(-1 + Power(x,2))*
              Cos(Pi*x) + x*Sin(Pi*x)
              )*Sin(4*Pi*y) + 
              45*Pi*Sin(Pi*t)*
              (x*
              Power(Sin(2*Pi*x),2) + 
              Pi*(-1 + Power(x,2))*
              Sin(4*Pi*x))*
              Sin(6*Pi*y)) - 
           ((8*Power(Pi,2)*
              (-1 + Power(x,2))*
              Power(Cos(Pi*x),2)*
              Sin((Pi*t)/2.)*
              Power(Sin(2*Pi*y),2) - 
              8*Sin((Pi*t)/2.)*
              Sin(Pi*x)*
              (-4*Pi*x*Cos(Pi*x) + 
              (-1 + 
              Power(Pi,2)*
              (-1 + Power(x,2)))*
              Sin(Pi*x))*
              Power(Sin(2*Pi*y),2) + 
              15*Sin(Pi*t)*
              (1 + 
              (-1 + 
              8*Power(Pi,2)*
              (-1 + Power(x,2)))*
              Cos(4*Pi*x) + 
              8*Pi*x*Sin(4*Pi*x))*
              Power(Sin(3*Pi*y),2))*
              (8*Pi*
              (-1 + Power(x,2))*
              Sin((Pi*t)/2.)*
              Power(Sin(Pi*x),2)*
              Sin(4*Pi*y) + 
              5*
              (-40 + 
              9*Pi*(-1 + Power(x,2))*
              Sin(Pi*t)*
              Power(Sin(2*Pi*x),2)*
              Sin(6*Pi*y))))/4.))/
	  100.)/100.);
}





int main(int argc, char** argv)
{
  try
    {
      Time timer("total");
      timer.start();
      int nx = 16;
      Sundance::setOption("nx", nx, "grid size");

      int nSteps = 16;
      Sundance::setOption("nSteps", nSteps, "number of timesteps");

      double tCur = 0.0; // time can't hit 0 + pi*k since q has a csc(t) term
      Sundance::setOption("iInit", tCur, "initial time");
      
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
      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, nx, 1,
                                                         0.0, 1.0, nx, 1,
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

      RCP<ProblemDescription> prob = rcp(new ProblemDescription(mesh,Stokes));
      prob->addCondition(BodyForce, interior, StateDependentFunction(qExFunc));
      
      // Create no slip Boundary Conditions on top and bottom
      prob->addCondition(NoSlip,top,StateDependentFunction(List(1.0,0.0)));
      prob->addCondition(NoSlip,bottom,StateDependentFunction(List(1.0,0.0)));
      
      //Set inflow and outflow conditions
      prob->addCondition(NoSlip,left,StateDependentFunction(List(1.0,0.0)));
      prob->addCondition(Traction,right,StateDependentFunction(0.0));

      // Set IC
      prob->setIC(StateDependentFunction(List(1.0,0.0)),StateDependentFunction(pExFunc));
      prob->setNuEff(1.0); // What is NuEff? Kinematic viscosity needs to be 1
      

      PCDControlParams pcdParams;
      NewtonControlParams newtonParams;

      SteppingType method = FullyImplicit;
      //SteppingType method = SemiImplicit;
      //pcdParams.FParams.method = "CG";

      PCDStepper stepper(prob, pcdParams, newtonParams, order, method);
      RCP<FlowState> state = stepper.state();

      double dt = (tFinal-tCur)/nSteps;
      
      state->initialize(tCur, dt, uExFunc, pExFunc);

      string outDir = "Results/Stokes/";
	    //string outDir = "Results/";
   //   system( ("mkdir -p " + outDir).c_str() ); Done in DefaultOutputManager constructor
      string filename = outDir + "forward_problem_TransientChannel_nx" + Teuchos::toString(nx)
	+ "-nt-" + Teuchos::toString(nSteps);

      system( ("rm -fr " + filename).c_str() );

      ErrorChecker check(state, uExFunc, pExFunc);

//Added by me to write out the snapshots
      RCP<DefaultOutputManager> output
	= rcp(new DefaultOutputManager(filename, "st", new ExodusWriterFactory()));

      // Same logic as VientoUniformStepController.cpp; I have taken it
      //out so I can utilize ErrorChecker
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
	  FieldWriter w = new VTKWriter(filename + "/" + "step"+ Teuchos::toString(i));
	  check.write(w);
	  output->write(i, stepper.state(), 1);
	}

      timer.stop();
      Out::root() << "nx=" << nx << " nt=" << nSteps << " uErr=" << check.uErrL2()
		  << " pErr=" << check.pErrL2() << " runtime=" << timer.totalElapsedTime()
		  << endl << endl << endl;

    }
  catch(std::exception& e)
    {
      Sundance::handleException(e);
    }
  Sundance::finalize(); return Sundance::testStatus(); 
}
