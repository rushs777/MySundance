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

//For using newton-armijo
#include "PlayaDenseLUSolver.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "PlayaNewtonArmijoSolverImpl.hpp"
#include "PlayaNonlinearSolver.hpp"
#include "PlayaNonlinearSolverBuilder.hpp"


//Local files
#include "PlayaSVD.hpp"
#include "MMSQuadODE.hpp"
#include "MyNLO.hpp"
#include "MathematicaConverter.hpp"
#include "PMB.hpp"
#include "velocityROM.hpp"

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

      int nSteps = 16;
      Sundance::setOption("nSteps", nSteps, "Number of time steps");

      double tol = .999;
      Sundance::setOption("tol", tol, "Tolerance requirement for the number of basis functions to keep");

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

      VectorType<double> time_vecType = new SerialVectorType();
      VectorSpace<double> time_vecSpace = time_vecType.createEvenlyPartitionedSpace(MPIComm::self(), nSteps+1.0);

      //Needed for the integral
      CellFilter interior = new MaximalCellFilter();
      QuadratureFamily quad = new GaussianQuadrature(6);

      /***********************Build uRO **********************/
      
      // Read the snapshots into a matrix
      string outDir = "Results/forward_problem_TransientChannel_nx" + Teuchos::toString(nx) + "-nt-" + Teuchos::toString(nSteps);
      string velocityTag = "st-v";
      string filename = outDir + "/" + velocityTag;


      // Create a BasisFamily to express our solution
      Array<Sundance::BasisFamily> ubasis = List(new Sundance::Lagrange(2), new Sundance::Lagrange(2)); // 2nd order Piece-Wise Quad Lagrange in 2D
      // Define our vector type
      Playa::VectorType<double> epetraVecType = new Playa::EpetraVectorType();
      Sundance::DiscreteSpace velocityDS(mesh, ubasis, epetraVecType);

      string NLParamFile = "playa-newton-armijo.xml";     

      velocityROM ROM(filename, NLParamFile, velocityDS, u0, q, t, nSteps, deltat, tol, verbosity);
      ROM.initialize();
      ROM.generate_alpha();
      Array<Expr> uRO(ROM.get_uRO() );
      //Array<Expr> phi(ROM.get_phi() ); // These are the POD basis functions

      Vector<double> uError = time_vecSpace.createMember();
      for(int time = 0; time < uRO.length(); time++)
	{
	  t.setParameterValue(tInit+time*deltat);
	  uError[time] = L2Norm(mesh, interior, uExact - uRO[time], quad);
	  SUNDANCE_ROOT_MSG2(verbosity, "Error for uRO at time " + Teuchos::toString(tInit+time*deltat) + "= " + Teuchos::toString(uError[time]));
	}
      /***********************end of uRO **********************/


      

      /*******************Begin trying to do PMB****************************************/

      // BasisFamily for p
      BasisFamily pbasis = new Sundance::Lagrange(1);
      DiscreteSpace pressureDS(mesh, pbasis, epetraVecType);

      // mesh.spatialDim() returns n for nD
      int dim = mesh.spatialDim();

      // Define our differential operators; note Derivative(x=0)
      Expr grad = gradient(dim);


      /*******************Generate pRO****************************************/

      string pressureTag = "st-p";
      string pressureFilename = outDir + "/" + pressureTag;

      PMB pmb(pressureFilename, uRO, q, t, pressureDS, deltat, tol, verbosity);
      pmb.initialize();
      pmb.generate_beta();
      Array<Expr> psi(pmb.get_psi() );

      SUNDANCE_ROOT_MSG1(verbosity, "Staring to build reduced-order p from beta");
      Array<Expr> pRO(pmb.get_pRO() );

      Vector<double> pError = time_vecSpace.createMember();

      for(int time = 0; time < pRO.length(); time++)
	{
	  t.setParameterValue(tInit + time*deltat);
	  pError[time] = L2Norm(mesh, interior, pExact - pRO[time], quad);
	  SUNDANCE_ROOT_MSG2(verbosity, "Error for pRO at time " + Teuchos::toString(tInit+time*deltat) + "= " + Teuchos::toString(pError[time]));
	}
      /***********************End of pRO **********************/
	    

      /***********************Build beta through projection **********************/
      
      //Find the exact betas
      Array<Vector<double> > beta(nSteps+1);

      // Based off the value for Rp, create an appropriate VectorSpace<double>
      int Rp = psi.length();
      VectorType<double> Rp_vecType = new SerialVectorType();
      VectorSpace<double> Rp_vecSpace = Rp_vecType.createEvenlyPartitionedSpace(MPIComm::self(), Rp);
      for(int count = 0; count<beta.length(); count++)
	beta[count] = Rp_vecSpace.createMember();

      // Create pbar
      LinearOperator<double> Qprime = snapshotToMatrix(pressureFilename, nSteps, mesh);

      // Calculate pbar(x)
      Vector<double> pbarVec = Qprime.range().createMember();
      Vector<double> ones = Qprime.domain().createMember();
      ones.setToConstant(1.0);
      Qprime.apply(ones, pbarVec);
      pbarVec *= (1.0/ (nSteps+1.0) );
      Expr pbar = new DiscreteFunction(pressureDS, serialToEpetra(pbarVec));

      for(int tIndex=0; tIndex < beta.length(); tIndex++)
	{
	  t.setParameterValue(tInit + tIndex*deltat);
	  for(int r=0; r<R; r++)
	    {
	      // beta_r(t_m) = ( pEx(t_m, x, y) - pbar(x,y), psi[r] )
	      FunctionalEvaluator ExactEvaluator = FunctionalEvaluator(mesh, Integral(interior, (pExact-pbar)*psi[r], quad));
	      beta[tIndex][r] = ExactEvaluator.evaluate();
	    }

	}

      // Build pIsh
      Array<Expr> pIsh(nSteps+1, pbar);
      for(int time=0; time < (nSteps+1); time++)
	{
	  for(int r = 0; r < R; r++)
	    pIsh[time] = pIsh[time] + beta[time][r]*psi[r];
	}

      SUNDANCE_ROOT_MSG1(verbosity, "Comparing pExact(t_n) to pIsh(t_n)");
      Vector<double> errorIsh = time_vecSpace.createMember();
            
      for(int time = 0; time < pRO.length(); time++)
	{
	  t.setParameterValue(tInit + time*deltat);
	  errorIsh[time] = L2Norm(mesh, interior, pExact - pIsh[time], quad);
	  SUNDANCE_ROOT_MSG2(verbosity, "Error for pIsh at time " + Teuchos::toString(tInit+time*deltat) + "= " + Teuchos::toString(errorIsh[time]));
	}

      cout << "||pExact - pIsh||_2:\t "  << errorIsh.norm2() << endl;
      cout << "||pExact - pIsh||_inf:\t " << errorIsh.normInf() << endl;

      /***********************End of build beta through projection **********************/
      

 
      
      // Visualize the results
      SUNDANCE_ROOT_MSG1(verbosity, "Writing results to file");
      string vtkDir = "Results/ROM/pRO/";
      string vtkfilename = "nx"+Teuchos::toString(nx)+"nt"+Teuchos::toString(nSteps);
      vtkDir = vtkDir + vtkfilename + "/";
      system( ("mkdir -p " + vtkDir).c_str() ); 

      DiscreteSpace scalarDS(mesh, new Lagrange(1), new EpetraVectorType());
      L2Projector uEx_projector(velocityDS, uExact);
      L2Projector pEx_projector(pressureDS, pExact);
      //Write uExact for all the time steps
      for(int time=0; time< nSteps+1; time++)
	{
	  FieldWriter writer = new VTKWriter(vtkDir+vtkfilename+"step"+Teuchos::toString(time));
	  writer.addMesh(mesh);
	  t.setParameterValue(tInit+time*deltat);
	  L2Projector uRO_projector(velocityDS, uRO[time]);
	  L2Projector uErrorProjector(velocityDS, uExact - uRO[time]);
	  L2Projector pRO_projector(pressureDS, pRO[time]);
	  L2Projector pErrorProjector(pressureDS, pExact - pRO[time]);

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
	  writer.addField("uExact[0]", new ExprFieldWrapper(uEx_projector.project()[0]) );
	  writer.addField("uExact[1]", new ExprFieldWrapper(uEx_projector.project()[1]) );
	  writer.addField("pExact", new ExprFieldWrapper(pEx_projector.project()) );
	  writer.addField("uRO[0]", new ExprFieldWrapper(uRO_projector.project()[0]) );
	  writer.addField("uRO[1]", new ExprFieldWrapper(uRO_projector.project()[1]) );
	  writer.addField("pRO", new ExprFieldWrapper(pRO_projector.project()) );
	  writer.addField("uError[0]", new ExprFieldWrapper(uErrorProjector.project()[0]) );
	  writer.addField("uError[1]", new ExprFieldWrapper(uErrorProjector.project()[1]) );
	  writer.addField("pError", new ExprFieldWrapper(pErrorProjector.project()) );
	  writer.write();
	}

      /*
      CellFilter bdry = new BoundaryCellFilter();
      DiscreteSpace bdryDS(mesh, new Lagrange(1), bdry, new EpetraVectorType());
      string vtkCheck = "ngradp";
      Expr nHat = CellNormalExpr(mesh.spatialDim(), "nHat");	  
      FieldWriter pressureWriter = new VTKWriter(vtkDir+vtkCheck);
      pressureWriter.addMesh(mesh);
      
      Expr ngradpExact = nHat*(grad*pExact);
      cout << "Here is ngradp " << ngradpExact << endl;
      L2Projector ngradpExact_projector(bdryDS, ngradpExact); // Code throws seg fault here
      pressureWriter.addField("ngradpExact", new ExprFieldWrapper(ngradpExact_projector.project()));

      cout << "hi " << endl;
      for(int r=0; r<R; r++)
	{
	  Expr ngradpsi = nHat*(grad*psi[r]);
	  L2Projector ngradpsi_projector(bdryDS, ngradpsi);
	  pressureWriter.addField("ngradpsi["+Teuchos::toString(r)+"]", new ExprFieldWrapper(ngradpsi_projector.project()));
      cout << "bye" << endl;
	}


      Expr ngradpBar = nHat*(grad*pbar);
      L2Projector ngradpBar_projector(bdryDS, ngradpBar);
      pressureWriter.addField("ngradpBar", new ExprFieldWrapper(ngradpBar_projector.project()));
      
      pressureWriter.write();
      */
      


      cout << "nt = " << nSteps << "\t nx = " << nx << endl;
      Out::root() << "||2norm of the error in u at all timesteps||_2  :\t " << uError.norm2() << endl;
      Out::root() << "||2norm of the error in u at all timesteps||_inf :\t " << uError.normInf() << endl;
      cout << "||2norm of the error in p at all timesteps||_2 :\t "  << pError.norm2() << endl;
      cout << "||2norm of the error in p at all timesteps||_inf :\t " << pError.normInf() << endl;
      cout << "Number of modes kept for pressure: " << Rp << endl;

      timer.stop();
      cout << "runtime=" << timer.totalElapsedTime() << endl;
    }
  catch(std::exception& e)
    {
      Out::root() << "Caught exception: " << e.what() << std::endl;
      return -1;
    }
}

