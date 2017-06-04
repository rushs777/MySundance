#include "Teuchos_GlobalMPISession.hpp"
#include "PlayaSerialVectorType.hpp"
#include "PlayaDenseSerialMatrix.hpp"
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaSerialEpetraAdapter.hpp" // Allows me to convert between Epetra and Serial Vectors
#include "PlayaOut.hpp"
#include <vector>


#include "VientoSnapshotIO.hpp" //For snapshotToMatrix
#include "denseSerialMatrixIO.hpp"
#include "QuadraticODERHSBase.hpp"
#include "MathematicaConverter.hpp"
#include "MMSQuadODE.hpp"
#include "MyNLO.hpp"
#include "velocityROM.hpp"

//For using newton-armijo
#include "PlayaDenseLUSolver.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "PlayaNewtonArmijoSolverImpl.hpp"
#include "PlayaNonlinearSolver.hpp"
#include "PlayaNonlinearSolverBuilder.hpp"


//Local files
#include "PlayaSVD.hpp"

#include "Sundance.hpp"

using std::vector;
using std::cout;
//using std::endl;
using namespace Teuchos;
using namespace Playa;
using namespace PlayaExprTemplates;




int main(int argc, char *argv[]) 
{
  try
    {
      Time timer("total");
      timer.start();
      
      int verbosity = 1;	
      Sundance::setOption("verbosity", verbosity, "verbosity level");
      
      int nx = 16;
      Sundance::setOption("nx", nx, "Number of elements along each axis");

      double tol = .999;
      Sundance::setOption("tol", tol, "Tolerance requirement for the number of basis functions to keep");      
      int nSteps = 16;
      Sundance::setOption("nSteps", nSteps, "Number of time steps");

      double tInit = 0.0; 
      Sundance::setOption("tInit", tInit, "initial time");

      double tFinal = 1.0;
      Sundance::setOption("tFinal", tFinal, "final time");

      bool AreMatrixAndTensorInFile = false;
      Sundance::setOption("MatrixAndTensorInFile", "MatrixAndTensorNotInFile", AreMatrixAndTensorInFile, "true if the matrix and tensor are available text files. false if they need to be created");

      Sundance::init(&argc, &argv);
      GlobalMPISession session(&argc, &argv);

      // Define our coordinates
      Expr x = new CoordExpr(0,"x");
      Expr y = new CoordExpr(1,"y");
      Expr t = new Sundance::Parameter(tInit);
      double deltat = (tFinal-tInit)/nSteps;
      double R = 1.0;

      // Define initial conditions, i.e. u(x,y,0)
      Expr u0 = List(1.0,0.0);

      // Define the forcing term
      Expr q = List((-(Power(Pi,2)*
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

      // Define our solution
      Expr uExact = List(1 - (Pi*(-1 + Power(x,2))*
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

      Expr pExact = 0.0;


      

      // Define our mesh
      MeshType meshType = new Sundance::BasicSimplicialMeshType();
      double xmin = 0.0;
      double xmax = 1.0;
      double ymin = 0.0;
      double ymax = 1.0;
      MeshSource mesher = new Sundance::PartitionedRectangleMesher(xmin,xmax,nx,1,ymin,ymax,nx,1,meshType,0);
      Mesh mesh = mesher.getMesh();



      // Read the snapshots into a matrix
      string outDir = "Results/forward_problem_TransientChannel_nx" + Teuchos::toString(nx) + "-nt-" + Teuchos::toString(nSteps);
      string tag = "st-v";
      string filename = outDir + "/" + tag;

      // Create a BasisFamily to express our solution
      Array<Sundance::BasisFamily> ubasis = List(new Sundance::Lagrange(2), new Sundance::Lagrange(2)); // 2nd order Piece-Wise Quad Lagrange in 2D
      // Define our vector type
      Playa::VectorType<double> epetraVecType = new Playa::EpetraVectorType();
      Sundance::DiscreteSpace velocityDS(mesh, ubasis, epetraVecType);

      string NLParamFile = "playa-newton-armijo.xml";

      // Create a velocityROM object
      velocityROM ROM(filename, NLParamFile, velocityDS, u0, q, t, nSteps, deltat, tol, verbosity);
      ROM.initialize();
      ROM.generate_alpha();
      Array<Expr> uRO(ROM.get_uRO() );
     
      VectorType<double> time_vecType = new SerialVectorType();
      VectorSpace<double> time_vecSpace = time_vecType.createEvenlyPartitionedSpace(MPIComm::self(), nSteps+1.0);

      //Needed for the integral
      CellFilter interior = new MaximalCellFilter();
      QuadratureFamily quad = new GaussianQuadrature(6);


      // Compute alpha "exactly"
      Array<Expr> phi(ROM.get_phi() ); // These are the POD basis functions
      



      // Based off the value for Ru, create an appropriate VectorSpace<double>
      int Ru = phi.length();
      VectorType<double> R_vecType = new SerialVectorType();
      VectorSpace<double> R_vecSpace = R_vecType.createEvenlyPartitionedSpace(MPIComm::self(), Ru);
      
      //Find the exact alphas
      Array<Vector<double> > alpha(nSteps+1);
      for(int count = 0; count<alpha.length(); count++)
	alpha[count] = R_vecSpace.createMember();



      // Calculate ubar(x)
      LinearOperator<double> Wprime = snapshotToMatrix(filename, nSteps, mesh);
      SUNDANCE_ROOT_MSG2(verbosity, "Size of W: " << Wprime.range().dim() << " by " << Wprime.domain().dim());
      Vector<double> ubarVec = Wprime.range().createMember();
      Vector<double> ones = Wprime.domain().createMember();
      ones.setToConstant(1.0);
      Wprime.apply(ones, ubarVec);
      ubarVec *= (1.0/ (nSteps+1.0) );
      Expr ubar = new DiscreteFunction(velocityDS, serialToEpetra(ubarVec));

      //Find the exact alphas
      for(int tIndex=0; tIndex < alpha.length(); tIndex++)
	{
	  t.setParameterValue(tInit+tIndex*deltat);
	  for(int r=0; r<Ru; r++)
	    {
	      // alpha_r(t_m) = ( uEx(t_m, x, y), phi[r] )
	      //FunctionalEvaluator ExactEvaluator(mesh, Integral(interior, (uExact-ubar)*phi[r], quad));
	      FunctionalEvaluator ExactEvaluator = FunctionalEvaluator(mesh, Integral(interior, (uExact-ubar)*phi[r], quad));
	      alpha[tIndex][r] = ExactEvaluator.evaluate();
	    }
	}

      Array<Vector<double> > soln(ROM.get_alpha() ); // Approximate alphas


      
      Vector<double> alphaError = time_vecSpace.createMember();
      for(int i=0; i < alpha.length(); i++)
	{
	  alphaError[i] = (alpha[i] - soln[i]).norm2();
	  if(verbosity>=2)
	    cout << "Error for alpha(t=" << i << "): " << alphaError[i] << endl;

	  if(verbosity>=3)
	    {
	      cout << "Exact alpha(t=" << i << "): " << endl << alpha[i] << endl;	
	      cout << "Approximate alpha(t=" << i << "): " << endl << soln[i] << endl << endl;
	    }
	}

      cout << "Run for nx = " << nx << ", nSteps = " << nSteps << endl;
      cout << "||2norm of the error in alpha at all timesteps||_2  :\t "  << alphaError.norm2() << endl;
      cout << "||2norm of the error in alpha at all timesteps||_inf:\t " << alphaError.normInf() << endl;
      

      SUNDANCE_ROOT_MSG2(verbosity, "Comparing uExact(t_n) to uRO(t_n)");
      Vector<double> uError = time_vecSpace.createMember();
      
      for(int time = 0; time < uRO.length(); time++)
	{
	  t.setParameterValue(tInit+time*deltat);
	  uError[time] = L2Norm(mesh, interior, uExact - uRO[time], quad);
	  double uNorm = L2Norm(mesh, interior, uExact, quad);
	  /* print the relative error */
	  SUNDANCE_ROOT_MSG2(verbosity, "Relative Error for uRO at time " << tInit+time*deltat << " = " << uError[time]/uNorm);
	}
      
      timer.stop();
  SUNDANCE_ROOT_MSG1(verbosity, "Number of velocity modes kept: " + Teuchos::toString(ROM.get_phi().size()));
      Out::root() << "||2norm of the error in u at all timesteps||_2  :\t " << uError.norm2() << endl;
      Out::root() << "||2norm of the error in u at all timesteps||_inf:\t " << uError.normInf() << endl;

     
      // Visualize the results
      SUNDANCE_ROOT_MSG1(verbosity, "Writing results to file");
      string vtkDir = "Results/ROM/uRO/";
      string vtkfilename = "nx"+Teuchos::toString(nx)+"nt"+Teuchos::toString(nSteps);
      vtkDir = vtkDir + vtkfilename + "/";
      system( ("mkdir -p " + vtkDir).c_str() ); 

      DiscreteSpace scalarDS(mesh, new Lagrange(1), new EpetraVectorType());
      L2Projector projector(velocityDS, uExact);
      //Write uExact for all the time steps
      for(int time=0; time< nSteps+1; time++)
	{
	  FieldWriter writer = new VTKWriter(vtkDir+vtkfilename+"step"+Teuchos::toString(time));
	  writer.addMesh(mesh);
	  t.setParameterValue(time*deltat);
	  L2Projector projectorRO(velocityDS, uRO[time]);
	  L2Projector uErrorProjector(velocityDS, uExact - uRO[time]);
	  Expr absErr = sqrt( (uExact - uRO[time])*(uExact - uRO[time]));
	  Expr absU = sqrt(uExact * uExact);
	  Expr absURO = sqrt(uRO[time] * uRO[time]);
	  L2Projector uMagProj(scalarDS, absU);
	  L2Projector uROMagProj(scalarDS, absURO);
	  L2Projector absErrorProj(scalarDS, absErr);
	  L2Projector relErrorProj(scalarDS, absErr / (absU + 1.0));
	  writer.addField("uMag", new ExprFieldWrapper(uMagProj.project()[0]) );
	  writer.addField("uROMag", new ExprFieldWrapper(uROMagProj.project()[0]) );
	  writer.addField("errAbs", new ExprFieldWrapper(absErrorProj.project()[0]) );
	  writer.addField("errRel", new ExprFieldWrapper(relErrorProj.project()[0]) );
	  writer.addField("uExact[0]", new ExprFieldWrapper(projector.project()[0]) );
	  writer.addField("uExact[1]", new ExprFieldWrapper(projector.project()[1]) );
	  writer.addField("uRO[0]", new ExprFieldWrapper(projectorRO.project()[0]) );
	  writer.addField("uRO[1]", new ExprFieldWrapper(projectorRO.project()[1]) );
	  writer.addField("uError[0]", new ExprFieldWrapper(uErrorProjector.project()[0]) );
	  writer.addField("uError[1]", new ExprFieldWrapper(uErrorProjector.project()[1]) );
      	  writer.write();	  
	}

      timer.stop();
      Out::root() << "runtime=" << timer.totalElapsedTime() << endl << endl;
	
    }
  catch(std::exception& e)
    {
      Out::root() << "Caught exception: " << e.what() << std::endl;
      return -1;
    }
}

