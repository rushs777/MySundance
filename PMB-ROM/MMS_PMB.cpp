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
      string outDir = "/home/sirush/PhDResearch/ODETest/Results";
      string fileDir = outDir + "/pcd-err-nx-" + Teuchos::toString(nx)
	+ "-nt-" + Teuchos::toString(nSteps);
      string tag = "st-v";
      string filename = fileDir + "/" + tag;
      Playa::LinearOperator<double> W = snapshotToMatrix(filename, nSteps, mesh);
      SUNDANCE_ROOT_MSG2(verbosity, "Size of W: " << W.range().dim() << " by " << W.domain().dim());

      // Perform the POD of the matrix W
      Playa::LinearOperator<double> U;
      Playa::LinearOperator<double> Phi;
      Playa::Vector<double> sigma;

      // Create a BasisFamily to express our solution
      Array<Sundance::BasisFamily> ubasis = List(new Sundance::Lagrange(2), new Sundance::Lagrange(2)); // 2nd order Piece-Wise Quad Lagrange in 2D
      // Define our vector type
      Playa::VectorType<double> epetraVecType = new Playa::EpetraVectorType();
      Sundance::DiscreteSpace ds(mesh, ubasis, epetraVecType);

      // W and mesh need to be defined
      SUNDANCE_ROOT_MSG1(verbosity, "Entering POD");
      POD(W,sigma,U,Phi,ds,verbosity);
      SUNDANCE_ROOT_MSG2(verbosity, "POD finished");

      Vector<double> lambda = sigma.copy();
      for(int count = 0; count < sigma.dim(); count++)
	lambda[count] = sqrt(lambda[count]);
      double lambdaTotal = lambda.norm1();
      double lambdaSum = 0.0;
      int R = 0;
      for(int i = 0; i<lambda.dim(); i++)
	{
	  lambdaSum += lambda[i];
	  if(lambdaSum/lambdaTotal >= .999)
	    {
	      // R is the number of lambdas to keep
	      R = i + 1;
	      SUNDANCE_ROOT_MSG2(verbosity, "Number of lambda[i] kept: " + Teuchos::toString(R));
	      break;
	    }
	}

      // Based off the value for R, create an appropriate VectorSpace<double>
      VectorType<double> R_vecType = new SerialVectorType();
      VectorSpace<double> R_vecSpace = R_vecType.createEvenlyPartitionedSpace(MPIComm::self(), R);




      // Looking at PlayaSVD.cpp, Phi is a DenseSerialMatrix
      Playa::Vector<double> ej = Phi.domain().createMember();
      Playa::Vector<double> phiCoeff = Phi.range().createMember(); // These are the coefficient vectors
      Array<Expr> phi(R); // These are the POD basis functions
      SUNDANCE_ROOT_MSG2(verbosity, "Size of phi: " + Teuchos::toString(phi.size()));

      //Needed for the integral
      CellFilter interior = new MaximalCellFilter();
      QuadratureFamily quad4 = new GaussianQuadrature(4);

      // Get the Expr phi_r(x)
      for(int r = 0; r<R; r++)
	{
	  ej.zero();
	  ej[r] = 1.0;
	  phiCoeff.zero();
	  Phi.apply(ej,phiCoeff);
	  phi[r] = new DiscreteFunction(ds, serialToEpetra(phiCoeff)); //DiscreteFunction requires Epetra vectors
	  TEUCHOS_TEST_FOR_EXCEPTION( fabs(L2Norm(mesh, interior, phi[r], quad4) - 1.0) >= 1.0e-6,
				     runtime_error, "||phi["+Teuchos::toString(r)+"]|| = " + Teuchos::toString(L2Norm(mesh, interior, phi[r], quad4)) + " != 1");
	}

      //Find the exact alphas
      Array<Vector<double> > alpha(nSteps+1);
      for(int count = 0; count<alpha.length(); count++)
	alpha[count] = R_vecSpace.createMember();

      // Needs to be of size nSteps+1;
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
      

      SUNDANCE_ROOT_MSG1(verbosity, "Staring to build reduced-order u from alphaExact ");
      Array<Expr> uRO(nSteps+1);
      for(int n=0; n<alpha.length(); n++)
	{	
	  uRO[n] = List(0.0, 0.0);
	  for(int r=0; r<R; r++)
	    {
	      uRO[n] = uRO[n] + alpha[n][r]*phi[r]; 
	    }
	}

      SUNDANCE_ROOT_MSG1(verbosity, "Comparing uExact(t_n) to uRO(t_n)");
      Vector<double> l2norm = Phi.domain().createMember();
      for(int time = 0; time < uRO.length(); time++)
	{
	  t.setParameterValue(time*deltat);
	  l2norm[time] = L2Norm(mesh, interior, uExact - uRO[time], quad4);
	  SUNDANCE_ROOT_MSG2(verbosity, "Error for uROExact at time " + Teuchos::toString(time*deltat) + "= " + Teuchos::toString(l2norm[time]));
	}

      double nu = 1.0;
      Expr q = List((-36*Cos(y)*Sin(t)*Sin(x) + (((44 + 30*Cos(t) + 8*Cos(4*t) + 30*Cos(5*t))*
						   Cos(x) + 9*((3 + 2*Cos(t) + 2*Cos(5*t))*Cos(3*x) + Cos(5*x)))*Power(Sin(x),3) + 9*Cos(6*t)*Cos(2*x)*Power(Sin(2*x),3))*Power(Sin(y),2) - 3*(8*nu*Cos(2*t)*(-1 + 2*Cos(2*x)) + 6*nu*Cos(3*t)*(-1 + 5*Cos(4*x)) + 8*Sin(2*t)*Power(Sin(x),2) + 9*Sin(3*t)*Power(Sin(2*x),2))*Sin(2*y))/36., (6*nu*(2*Cos(2*t)*(-1 + 2*Cos(2*y))*Sin(2*x) + 3*Cos(3*t)*(-4 + 5*Cos(2*y))*Sin(4*x)) - 18*Cos(x)*Sin(t)*Sin(y) + 3*(4*Sin(2*t)*Sin(2*x) + 9*Sin(3*t)*Sin(4*x))*Power(Sin(y),2) + 2*Cos(y)*(Cos(2*t)*(4*Cos(2*t) - 3*Cos(3*t)*(-3 - 6*Cos(2*x) + Cos(4*x)))*
Power(Sin(x),2) + 9*Power(Cos(3*t),2)*Power(Sin(2*x),2))*Power(Sin(y),3))/18.);
      
      // BasisFamily for p
      //BasisArray pbasis;
      //pbasis.push_back(new Sundance::Lagrange(1));
      BasisFamily pbasis = new Sundance::Lagrange(1);
      DiscreteSpace p_ds(mesh, pbasis, epetraVecType);

      // mesh.spatialDim() returns n for nD
      int dim = mesh.spatialDim();

      // Define our differential operators; note Derivative(x=0)
      Expr grad = gradient(dim);

      /*******************Begin trying to do PMB****************************************/

      string pressureTag = "st-p";
      string pressureFilename = fileDir + "/" + pressureTag;

      PMB pmb(pressureFilename, uRO, p_ds, deltat, .999, verbosity);
      pmb.initialize();
      pmb.generate_beta();

      //Find the exact betas
      /*      Array<Vector<double> > beta(nSteps+1);
      for(int count = 0; count<beta.length(); count++)
	beta[count] = pressure_R_vecSpace.createMember();

      t.setParameterValue(0.0);
      for(int tIndex=0; tIndex < beta.length(); tIndex++)
	{
	  for(int r=0; r<pressure_R; r++)
	    {
	      // beta_r(t_m) = ( pEx(t_m, x, y), psi[r] )
	      FunctionalEvaluator ExactEvaluator(mesh, Integral(interior, pExact*psi[r], quad4));
	      beta[tIndex][r] = ExactEvaluator.evaluate();
	    }
	  t.setParameterValue(t.getParameterValue()+deltat);
	}
      */






      SUNDANCE_ROOT_MSG1(verbosity, "Staring to build reduced-order p from beta");
      Array<Expr> pRO(pmb.get_pRO() );

      SUNDANCE_ROOT_MSG1(verbosity, "Comparing pExact(t_n) to pRO(t_n)");
      // This needs be fixed. We need a VectorSpace for the number of time nodes
      Vector<double> pressure_l2norm = Phi.domain().createMember();
      for(int time = 0; time < pRO.length(); time++)
	{
	  t.setParameterValue(time*deltat);
	  pressure_l2norm[time] = L2Norm(mesh, interior, pExact - pRO[time], quad4);
	  SUNDANCE_ROOT_MSG2(verbosity, "Error for pROPMB at time " + Teuchos::toString(time*deltat) + "= " + Teuchos::toString(pressure_l2norm[time]));
	}
 
      
      // Visualize the results
      string vtkDir = "Results/Visuals/pRO";
      system( ("mkdir -p " + vtkDir).c_str() ); 
      string vtkfilename = "nx"+Teuchos::toString(nx)+"nt"+Teuchos::toString(nSteps);
      FieldWriter writer = new VTKWriter(vtkDir+vtkfilename);
      writer.addMesh(mesh);

      L2Projector uEx_projector(ds, uExact);
      L2Projector pEx_projector(p_ds, pExact);
      
      //Write uExact for all the time steps
      for(int time=0; time<1; time++)
	{
	  t.setParameterValue(time*deltat);
	  L2Projector uRO_projector(ds, uRO[time]);
	  L2Projector pRO_projector(p_ds, pRO[time]);
	  L2Projector uErrorProjector(ds, uExact - uRO[time]);
	  L2Projector pErrorProjector(p_ds, pExact - pRO[time]);
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

      cout << "nt = " << nSteps << "\t nx = " << nx << endl;
      cout << "The 2-norm for the velocity error:\t " << l2norm.norm2() << endl;
      cout << "The 2-norm for the pROPMB error:\t "  << pressure_l2norm.norm2() << endl << endl;

	
    }
  catch(std::exception& e)
    {
      Out::root() << "Caught exception: " << e.what() << std::endl;
      return -1;
    }
}

