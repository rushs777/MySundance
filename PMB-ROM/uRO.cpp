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


#include <sys/resource.h>
void memcheck()
{
  Tabs tab;
  struct rusage r_usage;
  getrusage(RUSAGE_SELF, &r_usage);
  Out::root() << tab << "Memory used (MB): " << r_usage.ru_maxrss / ((double) 1.0e6) << endl;
}


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
      int verbosity = 1;	
      Sundance::setOption("verbosity", verbosity, "verbosity level");
      
      int nx = 64;
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

      // Define our solution
      Expr uExact = List((2*Cos(2*t)*Cos(y)*Power(Sin(x),2)*Sin(y))/3. + 
			 (Cos(3*t)*Cos(y)*Power(Sin(2*x),2)*Sin(y))/2.,
			 (-2*Cos(2*t)*Cos(x)*Sin(x)*Power(Sin(y),2))/3. - 
			 Cos(3*t)*Cos(2*x)*Sin(2*x)*Power(Sin(y),2));

      Expr pExact = Cos(x)*Cos(y)*Sin(t);

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
      string outDir = "../ODETest/Results";
      string fileDir = outDir + "/pcd-err-nx-" + Teuchos::toString(nx)
	+ "-nt-" + Teuchos::toString(nSteps);
      string tag = "st-v";
      string filename = fileDir + "/" + tag;
      /*
      Playa::LinearOperator<double> W = snapshotToMatrix(filename, nSteps, mesh);
      SUNDANCE_ROOT_MSG2(verbosity, "Size of W: " << W.range().dim() << " by " << W.domain().dim());

      // Perform the POD of the matrix W
      Playa::LinearOperator<double> U;
      Playa::LinearOperator<double> Phi;
      Playa::Vector<double> sigma;
      */

      // Create a BasisFamily to express our solution
      Array<Sundance::BasisFamily> ubasis = List(new Sundance::Lagrange(2), new Sundance::Lagrange(2)); // 2nd order Piece-Wise Quad Lagrange in 2D
      // Define our vector type
      Playa::VectorType<double> epetraVecType = new Playa::EpetraVectorType();
      Sundance::DiscreteSpace ds(mesh, ubasis, epetraVecType);

      string NLParamFile = "playa-newton-armijo.xml";

      // Based off the value for R, create an appropriate VectorSpace<double>
      int R = 2;
      VectorType<double> R_vecType = new SerialVectorType();
      VectorSpace<double> R_vecSpace = R_vecType.createEvenlyPartitionedSpace(MPIComm::self(), R);
      
      //Find the exact alphas
      Array<Vector<double> > alpha(nSteps+1);
      for(int count = 0; count<alpha.length(); count++)
	alpha[count] = R_vecSpace.createMember();

      alpha[0][0] = -2.18212;
      alpha[0][1] = -0.137505;      

      // Attempt to declare a velocityRO object
      velocityROM ROM(filename, NLParamFile, ds, alpha[0], nSteps, deltat, .999, verbosity);
      ROM.initialize();
      cout << "Value of alpha[0]" << endl << alpha[0] << endl;
      ROM.generate_alpha();
      Array<Expr> uRO(ROM.get_uRO() );
      Array<Expr> phi(ROM.get_phi() ); // These are the POD basis functions



      //Needed for the integral
      CellFilter interior = new MaximalCellFilter();
      QuadratureFamily quad4 = new GaussianQuadrature(4);

      //Find the exact alphas
      for(int tIndex=0; tIndex < alpha.length(); tIndex++)
	{
	  for(int r=0; r<R; r++)
	    {
	      // alpha_r(t_m) = ( uEx(t_m, x, y), phi[r] )
	      FunctionalEvaluator ExactEvaluator(mesh, Integral(interior, uExact*phi[r], quad4));
	      alpha[tIndex][r] = ExactEvaluator.evaluate();
	    }
	  t.setParameterValue(t.getParameterValue()+deltat);
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
      
      cout << "||alphaExact - alphaApprox||_2:\t\t "  << alphaError.norm2() << endl;
      cout << "||alphaExact - alphaApprox||_inf:\t\t " << alphaError.normInf() << endl;


      SUNDANCE_ROOT_MSG2(verbosity, "Comparing uExact(t_n) to uRO(t_n)");
      Vector<double> l2norm = time_vecSpace.createMember();
      for(int time = 0; time < uRO.length(); time++)
	{
	  t.setParameterValue(time*deltat);
	  l2norm[time] = L2Norm(mesh, interior, uExact - uRO[time], quad4);
	  SUNDANCE_ROOT_MSG2(verbosity, "Error for uRO at time " + Teuchos::toString(time*deltat) + "= " + Teuchos::toString(l2norm[time]));
	}
      Out::root() << "||uExact - uRO||_2:\t " << l2norm.norm2() << endl;
      Out::root() << "||uExact - uRO||_inf:\t " << l2norm.normInf() << endl << endl;

      SUNDANCE_ROOT_MSG1(verbosity, "Number of velocity modes kept: " + Teuchos::toString(phi.size()));

      cout << "Value for alpha(0) " << endl << soln[0] << endl;


      
      /*
      cout << "Staring to build reduced-order u from alphaExact " << endl;
      Array<Expr> uRO(nSteps+1);
      for(int n=0; n<alpha.length(); n++)
	{	
	  uRO[n] = List(0.0, 0.0);
	  for(int r=0; r<R; r++)
	    {
	      uRO[n] = uRO[n] + alpha[n][r]*phi[r]; 
	    }
	}

      cout << "Comparing uExact(t_n) to uRO(t_n)" << endl;
      Vector<double> l2norm = Phi.domain().createMember();
      cout << "size of l2norm: " << l2norm.dim() << endl;
      for(int time = 0; time < uRO.length(); time++)
	{
	  t.setParameterValue(time*deltat);
	  l2norm[time] = L2Norm(mesh, interior, uExact - uRO[time], quad4);
	  cout << "Error at time " << time*deltat << "= " << l2norm[time] << endl;
	}
      */

      
      
      // Visualize the results
      /*
      string vtkDir = "Results/Visuals/uRO/";
      string vtkfilename = "nx"+Teuchos::toString(nx)+"nt"+Teuchos::toString(nSteps);
      system( ("mkdir -p " + vtkDir).c_str() ); 
      FieldWriter writer = new VTKWriter(vtkDir+vtkfilename);
      writer.addMesh(mesh);

      L2Projector projector(ds, uExact);
      //Write uExact for all the time steps
      for(int time=0; time<1; time++)
	{
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
      */
	
    }
  catch(std::exception& e)
    {
      Out::root() << "Caught exception: " << e.what() << std::endl;
      return -1;
    }
}

