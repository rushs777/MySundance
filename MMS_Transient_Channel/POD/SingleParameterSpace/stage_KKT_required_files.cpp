#include "Sundance.hpp"
#include "VientoSnapshotIO.hpp"
#include "denseSerialMatrixIO.hpp"
#include "NSEProjectedODE.hpp"

#include "MathematicaConverter.hpp"

/*
 * The purpose of this file is to stage the necessary components before performing
 * ODE-constrained optimization by solving the KKT system. This program
 * creates the time-mean function, uB, and moves it along with the results of 
 * projecting the NSE onto the space spanned by phi, A, T, and b, into the POD
 * directory.
 */

int main(int argc, char* argv[])
{
  Time timer("total");
  timer.start();
  
  int nx = 8;
  Sundance::setOption("nx",nx,"Number of spatial elements along each axis");

  int nSteps = 8;
  Sundance:: setOption("nSteps",nSteps,"Number of time steps taken to get from tInit to tFinal");

  double tol = 0.999;
  Sundance::setOption("tol",tol,"tol is the tolerance used by the RIC to determine the number of basis functions to keep");

  int precision = 3;
  Sundance::setOption("precision",precision,"Number of significant digits to keep in a filename using tol"); 

  double Re = 1.0;
  Sundance::setOption("Re", Re, "Reynolds number");

  double U0 = 1.0;
  Sundance::setOption("U0", U0, "Characteristic velocity used for calculating nu");
  
  double L0 = 1.0;
  Sundance::setOption("L0", L0, "Characteristic length used for calculating nu");  
  
  int verbosity = 0;
  Sundance::setOption("verbosity",verbosity,"Level of verbosity for displayed output");

  double tFinal = 1.0;
  Sundance::setOption("tFinal",tFinal,"Final time value");

  Sundance::init(&argc, &argv);

  // Define our coordinates
  double tInit = 0.0;
  Expr x = new CoordExpr(0,"x");
  Expr y = new CoordExpr(1,"y");
  Expr t = new Sundance::Parameter(tInit);
  double R = Re; // Reynolds Number
  double deltat = (tFinal-tInit)/nSteps;

  // Define the forcing term
  Expr q =  List((Pi*(-((-1 + Power(x,2))*(32*Pi*Cos(4*Pi*y)*Sin((Pi*t)/2.)*Power(Sin(Pi*x),2) +270*Pi*Cos(6*Pi*y)*Sin(Pi*t)*Power(Sin(2*Pi*x),2))*(4*Sin((Pi*t)/2.)*Sin(Pi*x)*(Pi*(-1 + Power(x,2))*Cos(Pi*x) +x*Sin(Pi*x))*Power(Sin(2*Pi*y),2)+ 15*Sin(Pi*t)*Sin(2*Pi*x)*(2*Pi*(-1 + Power(x,2))*Cos(2*Pi*x) +x*Sin(2*Pi*x))*Power(Sin(3*Pi*y),2))) -100*Pi*(-1 + Power(x,2))*(4*Cos((Pi*t)/2.)*Power(Sin(Pi*x),2)*Sin(4*Pi*y) +45*Cos(Pi*t)*Power(Sin(2*Pi*x),2)*Sin(6*Pi*y)) +(200*(8*Power(Pi,2)*(-1 + Power(x,2))*Power(Cos(Pi*x),2)*Sin((Pi*t)/2.)*Sin(4*Pi*y) -8*Sin((Pi*t)/2.)*Sin(Pi*x)*(-4*Pi*x*Cos(Pi*x) +(-1 +9*Power(Pi,2)*(-1 + Power(x,2)))*Sin(Pi*x))*Sin(4*Pi*y) +(45*Sin(Pi*t)*(1 +18*Power(Pi,2) -18*Power(Pi,2)*Power(x,2) +(-1 +26*Power(Pi,2)*(-1 + Power(x,2)))*Cos(4*Pi*x) +8*Pi*x*Sin(4*Pi*x))*Sin(6*Pi*y))/2.))/R+ (8*Sin((Pi*t)/2.)*Sin(Pi*x)*(Pi*(-1 + Power(x,2))*Cos(Pi*x) +x*Sin(Pi*x))*Sin(4*Pi*y) +45*Sin(Pi*t)*Sin(2*Pi*x)*(2*Pi*(-1 + Power(x,2))*Cos(2*Pi*x) +x*Sin(2*Pi*x))*Sin(6*Pi*y))*(8*Pi*(-1 + Power(x,2))*Sin((Pi*t)/2.)*Power(Sin(Pi*x),2)*Sin(4*Pi*y) +5*(-40 +9*Pi*(-1 + Power(x,2))*Sin(Pi*t)*Power(Sin(2*Pi*x),2)*Sin(6*Pi*y)))))/20000.,(400*Pi*(2*Cos((Pi*t)/2.)*Sin(Pi*x)*(Pi*(-1 + Power(x,2))*Cos(Pi*x) +x*Sin(Pi*x))*Power(Sin(2*Pi*y),2) +15*Cos(Pi*t)*Sin(2*Pi*x)*(2*Pi*(-1 + Power(x,2))*Cos(2*Pi*x) +x*Sin(2*Pi*x))*Power(Sin(3*Pi*y),2)) -(800*Pi*(16*Pi*Power(Cos(2*Pi*y),2)*Sin((Pi*t)/2.)*Sin(Pi*x)*(Pi*(-1 + Power(x,2))*Cos(Pi*x) +x*Sin(Pi*x)) +135*Pi*Power(Cos(3*Pi*y),2)*Sin(Pi*t)*Sin(2*Pi*x)*(2*Pi*(-1 + Power(x,2))*Cos(2*Pi*x) +x*Sin(2*Pi*x)) +12*Pi*x*Power(Cos(Pi*x),2)*Sin((Pi*t)/2.)*Power(Sin(2*Pi*y),2)- 28*Pi*x*Sin((Pi*t)/2.)*Power(Sin(Pi*x),2)*Power(Sin(2*Pi*y),2)+ 6*Sin((Pi*t)/2.)*Sin(2*Pi*x)*Power(Sin(2*Pi*y),2)+ 12*Power(Pi,2)*Sin((Pi*t)/2.)*Sin(2*Pi*x)*Power(Sin(2*Pi*y),2)- 12*Power(Pi,2)*Power(x,2)*Sin((Pi*t)/2.)*Sin(2*Pi*x)*Power(Sin(2*Pi*y),2)+ 180*Pi*x*Power(Cos(2*Pi*x),2)*Sin(Pi*t)*Power(Sin(3*Pi*y),2)- 315*Pi*x*Sin(Pi*t)*Power(Sin(2*Pi*x),2)*Power(Sin(3*Pi*y),2)+ 45*Sin(Pi*t)*Sin(4*Pi*x)*Power(Sin(3*Pi*y),2)+ 255*Power(Pi,2)*Sin(Pi*t)*Sin(4*Pi*x)*Power(Sin(3*Pi*y),2)- 255*Power(Pi,2)*Power(x,2)*Sin(Pi*t)*Sin(4*Pi*x)*Power(Sin(3*Pi*y),2)))/R + 4*Pi*(4*Sin((Pi*t)/2.)*Sin(Pi*x)*(Pi*(-1 + Power(x,2))*Cos(Pi*x) +x*Sin(Pi*x))*Power(Sin(2*Pi*y),2) +15*Sin(Pi*t)*Sin(2*Pi*x)*(2*Pi*(-1 + Power(x,2))*Cos(2*Pi*x) +x*Sin(2*Pi*x))*Power(Sin(3*Pi*y),2))*(8*Sin((Pi*t)/2.)*Sin(Pi*x)*(Pi*(-1 + Power(x,2))*Cos(Pi*x) +x*Sin(Pi*x))*Sin(4*Pi*y) +45*Sin(Pi*t)*Sin(2*Pi*x)*(2*Pi*(-1 + Power(x,2))*Cos(2*Pi*x) +x*Sin(2*Pi*x))*Sin(6*Pi*y)) -(8*Power(Pi,2)*(-1 + Power(x,2))*Power(Cos(Pi*x),2)*Sin((Pi*t)/2.)*Power(Sin(2*Pi*y),2) -8*Sin((Pi*t)/2.)*Sin(Pi*x)*(-4*Pi*x*Cos(Pi*x) +(-1 +Power(Pi,2)*(-1 + Power(x,2)))*Sin(Pi*x))*Power(Sin(2*Pi*y),2) +15*Sin(Pi*t)*(1 +(-1 +8*Power(Pi,2)*(-1 + Power(x,2)))*Cos(4*Pi*x) +8*Pi*x*Sin(4*Pi*x))*Power(Sin(3*Pi*y),2))*(8*Pi*(-1 + Power(x,2))*Sin((Pi*t)/2.)*Power(Sin(Pi*x),2)*Sin(4*Pi*y) +5*(-40 +9*Pi*(-1 + Power(x,2))*Sin(Pi*t)*Power(Sin(2*Pi*x),2)*Sin(6*Pi*y))))/40000.);




  
  
  // Define our mesh
  // Might be a good idea to read this in from the ForwardProblem 
  MeshType spatialMeshType = new Sundance::BasicSimplicialMeshType();
  double xmin = 0.0;
  double xmax = 1.0;
  double ymin = 0.0;
  double ymax = 1.0;
  MeshSource spatialMesher = new Sundance::PartitionedRectangleMesher(xmin,xmax,nx,1,ymin,ymax,nx,1,spatialMeshType,0);
  Mesh spatialMesh = spatialMesher.getMesh();

  // Create a BasisFamily to express our solution; each spatial dim is Lagrange(2)
  Array<Sundance::BasisFamily> velBasis(spatialMesh.spatialDim());
  for(int i = 0; i < velBasis.length(); i++)
    velBasis[i] = new Sundance::Lagrange(2);
      
  // Define our vector type
  Playa::VectorType<double> epetraVecType = new Playa::EpetraVectorType();
  Sundance::DiscreteSpace velocityDS(spatialMesh, velBasis, epetraVecType);

  // Create the snapshot matrix W
  // Absolute Path
  //string outDir = "/home/sirush/PhDResearch/MMS_Transient_Channel/ForwardProblem/Results/Re";
  // Relative Path
  string outDir = "../../ForwardProblem/Results/Re";
  std::ostringstream ReynoldsString;
  ReynoldsString << std::setprecision(precision) << Re;
  outDir = outDir + ReynoldsString.str()
    + "/nx" + Teuchos::toString(nx) + "nt" + Teuchos::toString(nSteps) + "/";
  string snapshotLibraryFilePrefix = outDir + "st-v";
  LinearOperator<double> W = snapshotToMatrix(snapshotLibraryFilePrefix, nSteps, spatialMesh);
  
  // Create uB(x)
  Expr uB = timeMeanFunctionGenerator(W, velocityDS);

  // My current thinking is that this should go in the POD directory since it relies
  // on the same information
  string POD_DataDir = "Results/Re"
    + ReynoldsString.str() + "/nx" + Teuchos::toString(nx) + "nt" + Teuchos::toString(nSteps)
    + "/tol";
  std::ostringstream tolFileValue;
  tolFileValue << std::setprecision(precision) << tol;
  POD_DataDir = POD_DataDir + tolFileValue.str() + "/";
  string uB_storagePrefix = POD_DataDir + "uB";
  writeSnap(uB_storagePrefix,0,uB);

  // Also adding code that will move the creation of the RHS of projected NSE to here
  // These files will also be sent to the POD_DataDir
  // Read in the POD from file
  string POD_basis_fileprefix = POD_DataDir + "POD_basis";
  Array<Expr> phi;
  int Ru = 10000;
  for(int r = 0; r < Ru; r++)
    {
      try
	{
	  phi.push_back( readSnap(POD_basis_fileprefix, r, spatialMesh ) );
	}
      catch (std::runtime_error& e)
	{
	  Ru = r;
	}
    }

  // Calculate the kinematic viscosity nu
  double nu = (U0*L0)/Re;

  RCP<NSEProjectedODE> f = rcp(new NSEProjectedODE(phi, uB, q, t, deltat, nu, spatialMesh, false, verbosity));
  f->initialize();
  f->write_b(nSteps);

  // Move the files to the POD Directory
  int moveError = system( ("mv NSEProjectedODEFiles/* " + POD_DataDir + ".").c_str() );
  TEUCHOS_TEST_FOR_EXCEPTION( moveError == -1,
			      runtime_error,
			      "Failed to move A, T, or b to the POD Directory");


  timer.stop();
  cout << "runtime = " << timer.totalElapsedTime() << endl << endl << endl;


}
