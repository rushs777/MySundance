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
using std::endl;
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

      double tol = .999;
      Sundance::setOption("tol", tol, "Tolerance requirement for the number of basis functions to keep");
      
      int nx = 32;
      Sundance::setOption("nx", nx, "Number of elements along each axis");

      int nSteps = 32;
      Sundance::setOption("nSteps", nSteps, "Number of time steps");

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



      // Read the snapshots into a matrix
      //      string outDir = "/home/sirush/PhDResearch/ODETest/Results";
      string outDir = "Results/ForwardProblem/HG_Dirichlet_nx";
      string fileDir = outDir + Teuchos::toString(nx) + "_nt" + Teuchos::toString(nSteps);
      string tag = "st-v";
      string filename = fileDir + "/" + tag;

      // Create a BasisFamily to express our solution
      Array<Sundance::BasisFamily> ubasis = List(new Sundance::Lagrange(2), new Sundance::Lagrange(2)); // 2nd order Piece-Wise Quad Lagrange in 2D
      // Define our vector type
      Playa::VectorType<double> epetraVecType = new Playa::EpetraVectorType();
      Sundance::DiscreteSpace ds(mesh, ubasis, epetraVecType);

      string NLParamFile = "playa-newton-armijo.xml";     

      // Attempt to declare a velocityRO object
      velocityROM ROM(filename, NLParamFile, ds, u0, q, t, nSteps, deltat, tol, verbosity);
      ROM.initialize();
      ROM.generate_alpha();
      Array<Expr> uRO(ROM.get_uRO() );
      Array<Expr> phi(ROM.get_phi() ); // These are the POD basis functions



      //Needed for the integral
      CellFilter interior = new MaximalCellFilter();
      QuadratureFamily quad4 = new GaussianQuadrature(6);

      // Based off the value for R, create an appropriate VectorSpace<double>
      int R = phi.length();
      VectorType<double> R_vecType = new SerialVectorType();
      VectorSpace<double> R_vecSpace = R_vecType.createEvenlyPartitionedSpace(MPIComm::self(), R);
      
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
      Expr ubar = new DiscreteFunction(ds, serialToEpetra(ubarVec));

      //Find the exact alphas
      for(int tIndex=0; tIndex < alpha.length(); tIndex++)
	{
	  t.setParameterValue(0.0+tIndex*deltat);
	  for(int r=0; r<R; r++)
	    {
	      // alpha_r(t_m) = ( uEx(t_m, x, y) - ubar(x,y), phi[r] )
	      FunctionalEvaluator ExactEvaluator = FunctionalEvaluator(mesh, Integral(interior, (uExact-ubar)*phi[r], quad4));
	      alpha[tIndex][r] = ExactEvaluator.evaluate();
	    }
	}

      Array<Vector<double> > soln(ROM.get_alpha() );
      VectorType<double> time_vecType = new SerialVectorType();
      VectorSpace<double> time_vecSpace = time_vecType.createEvenlyPartitionedSpace(MPIComm::self(), nSteps+1.0);
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
      cout << "||alphaExact - alphaApprox||_2  :\t "  << alphaError.norm2() << endl;
      cout << "||alphaExact - alphaApprox||_inf:\t " << alphaError.normInf() << endl;


      SUNDANCE_ROOT_MSG2(verbosity, "Comparing uExact(t_n) to uRO(t_n)");
      Vector<double> l2norm = time_vecSpace.createMember();
      for(int time = 0; time < uRO.length(); time++)
	{
	  t.setParameterValue(time*deltat);
	  l2norm[time] = L2Norm(mesh, interior, uExact - uRO[time], quad4);
	  SUNDANCE_ROOT_MSG2(verbosity, "Error for uRO at time " + Teuchos::toString(time*deltat) + "= " + Teuchos::toString(l2norm[time]));
	}

      Out::root() << "||uExact - uRO||_2  :\t " << l2norm.norm2() << endl;
      Out::root() << "||uExact - uRO||_inf:\t " << l2norm.normInf() << endl;

      SUNDANCE_ROOT_MSG1(verbosity, "Number of velocity modes kept: " + Teuchos::toString(phi.size()));
      
     
      // Visualize the results
      SUNDANCE_ROOT_MSG1(verbosity, "Writing results to file");
      string vtkDir = "Results/ROM/uROFluctuations/";
      string vtkfilename = "nx"+Teuchos::toString(nx)+"nt"+Teuchos::toString(nSteps);
      vtkDir = vtkDir + vtkfilename + "/";
      system( ("mkdir -p " + vtkDir).c_str() ); 

      L2Projector projector(ds, uExact);
      //Write uExact for all the time steps
      for(int time=0; time< nSteps+1; time++)
	{
	  FieldWriter writer = new VTKWriter(vtkDir+vtkfilename+"step"+Teuchos::toString(time));
	  writer.addMesh(mesh);
	  t.setParameterValue(time*deltat);
	  L2Projector projectorRO(ds, uRO[time]);
	  L2Projector uErrorProjector(ds, uExact - uRO[time]);
	  writer.addField("uExact[0]", new ExprFieldWrapper(projector.project()[0]) );
	  writer.addField("uExact[1]", new ExprFieldWrapper(projector.project()[1]) );
	  writer.addField("uRO[0]", new ExprFieldWrapper(projectorRO.project()[0]) );
	  writer.addField("uRO[1]", new ExprFieldWrapper(projectorRO.project()[1]) );
	  writer.addField("uError[0]", new ExprFieldWrapper(uErrorProjector.project()[0]) );
	  writer.addField("uError[1]", new ExprFieldWrapper(uErrorProjector.project()[1]) );
      	  writer.write();	  
	}

	
    }
  catch(std::exception& e)
    {
      Out::root() << "Caught exception: " << e.what() << std::endl;
      return -1;
    }
}

