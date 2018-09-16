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

/*******************************************************
 *
 * NOTE: CODE IS HARDCODED FOR tInit TO BE 0
 *
 ******************************************************/


int main(int argc, char *argv[]) 
{
  try
    {
      Time timer("total");
      timer.start();
      
      int verbosity = 1;	
      Sundance::setOption("verbosity", verbosity, "verbosity level");
      
      int nx = 32;
      Sundance::setOption("nx", nx, "Number of elements along each axis");

      int nSteps = 32;
      Sundance::setOption("nSteps", nSteps, "Number of time steps");

      double tol = .999;
      Sundance::setOption("tol", tol, "Tolerance requirement for the number of basis functions to keep");

      double tFinal = 1.0;
      Sundance::setOption("tFinal", tFinal, "final time");

      bool AreMatrixAndTensorInFile = false;
      Sundance::setOption("MatrixAndTensorInFile", "MatrixAndTensorNotInFile", AreMatrixAndTensorInFile, "true if the matrix and tensor are available text files. false if they need to be created");

      Sundance::init(&argc, &argv);
      GlobalMPISession session(&argc, &argv);

      // Define our coordinates
      Expr x = new CoordExpr(0,"x");
      Expr y = new CoordExpr(1,"y");
      Expr t = new Sundance::Parameter(0.0);
      double deltat = tFinal/nSteps;
      double nu = 1.0;

      // Define our solution
      Expr uExact = List((2*Cos(2*t)*Cos(y)*Power(Sin(x),2)*Sin(y))/3. + 
			 (Cos(3*t)*Cos(y)*Power(Sin(2*x),2)*Sin(y))/2.,
			 (-2*Cos(2*t)*Cos(x)*Sin(x)*Power(Sin(y),2))/3. - 
			 Cos(3*t)*Cos(2*x)*Sin(2*x)*Power(Sin(y),2));

      Expr pExact = Cos(x)*Cos(y)*Sin(t);

      Expr u0 = List((2*Cos(y)*Power(Sin(x),2)*Sin(y))/3. + (Cos(y)*Power(Sin(2*x),2)*
                     Sin(y))/2., (-2*Cos(x)*Sin(x)*Power(Sin(y),2))/3. - Cos(2*x)*Sin(2*x)*
                     Power(Sin(y),2));

      Expr q = List((-36*Cos(y)*Sin(t)*Sin(x) + (((44 + 30*Cos(t) + 8*Cos(4*t) + 30*Cos(5*t))*
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

      // Define our mesh
      MeshType meshType = new Sundance::BasicSimplicialMeshType();
      double xmin = 0.0;
      double xmax = 2.0*Pi;
      double ymin = 0.0;
      double ymax = 2.0*Pi;
      MeshSource mesher = new Sundance::PartitionedRectangleMesher(xmin,xmax,nx,1,ymin,ymax,nx,1,meshType,0);
      Mesh mesh = mesher.getMesh();

      VectorType<double> time_vecType = new SerialVectorType();
      VectorSpace<double> time_vecSpace = time_vecType.createEvenlyPartitionedSpace(MPIComm::self(), nSteps+1.0);

      //Needed for the integral
      CellFilter interior = new MaximalCellFilter();
      QuadratureFamily quad = new GaussianQuadrature(6);

      /***********************Build uRO **********************/
      
      // Read the snapshots into a matrix
      string outDir = "Results/ForwardProblem/HG_Dirichlet_nx";
      string fileDir = outDir + Teuchos::toString(nx) + "_nt" + Teuchos::toString(nSteps);
      string velocityTag = "st-v";
      string filename = fileDir + "/" + velocityTag;


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

      Vector<double> l2norm = time_vecSpace.createMember();
      for(int time = 0; time < uRO.length(); time++)
	{
	  t.setParameterValue(time*deltat);
	  l2norm[time] = L2Norm(mesh, interior, uExact - uRO[time], quad);
	  SUNDANCE_ROOT_MSG2(verbosity, "Error for uRO at time " + Teuchos::toString(time*deltat) + "= " + Teuchos::toString(l2norm[time]));
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
      string pressureFilename = fileDir + "/" + pressureTag;

      PMB pmb(pressureFilename, uRO, q, t, pressureDS, deltat, tol, verbosity);
      pmb.initialize();
      pmb.generate_beta();
      Array<Expr> psi(pmb.get_psi() );

      SUNDANCE_ROOT_MSG1(verbosity, "Staring to build reduced-order p from beta");
      Array<Expr> pRO(pmb.get_pRO() );

      SUNDANCE_ROOT_MSG1(verbosity, "Comparing pExact(t_n) to pRO(t_n)");
      Vector<double> pressure_l2norm = time_vecSpace.createMember();

      for(int time = 0; time < pRO.length(); time++)
	{
	  t.setParameterValue(time*deltat);
	  pressure_l2norm[time] = L2Norm(mesh, interior, pExact - pRO[time], quad);
	  SUNDANCE_ROOT_MSG2(verbosity, "Error for pRO at time " + Teuchos::toString(time*deltat) + "= " + Teuchos::toString(pressure_l2norm[time]));
	}
      /***********************End of pRO **********************/
	    

      /***********************Build beta through projection **********************/
      
      //Find the exact betas
      Array<Vector<double> > beta(nSteps+1);

      // Based off the value for R, create an appropriate VectorSpace<double>
      int R = psi.length();
      VectorType<double> R_vecType = new SerialVectorType();
      VectorSpace<double> R_vecSpace = R_vecType.createEvenlyPartitionedSpace(MPIComm::self(), R);
      for(int count = 0; count<beta.length(); count++)
	beta[count] = R_vecSpace.createMember();

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
	  t.setParameterValue(0.0 + tIndex*deltat);
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

      Vector<double> errorIsh = time_vecSpace.createMember();
            
      for(int time = 0; time < pRO.length(); time++)
	{
	  t.setParameterValue(0.0 + time*deltat);
	  errorIsh[time] = L2Norm(mesh, interior, pExact - pIsh[time], quad);
	  SUNDANCE_ROOT_MSG2(verbosity, "Error for pIsh at time " + Teuchos::toString(time*deltat) + "= " + Teuchos::toString(errorIsh[time]));
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
	  t.setParameterValue(time*deltat);
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
      


      cout << "nt = " << nSteps << "\t nx = " << nx << endl;
      Out::root() << "||uExact - uRO||_2  :\t " << l2norm.norm2() << endl;
      Out::root() << "||uExact - uRO||_inf:\t " << l2norm.normInf() << endl;
      cout << "||pExact - pRO||_2:\t "  << pressure_l2norm.norm2() << endl;
      cout << "||pExact - pRO||_inf:\t " << pressure_l2norm.normInf() << endl;
      cout << "Number of modes kept for pressure: " << R << endl;

      timer.stop();
      cout << "runtime=" << timer.totalElapsedTime() << endl;
    }
  catch(std::exception& e)
    {
      Out::root() << "Caught exception: " << e.what() << std::endl;
      return -1;
    }
}

