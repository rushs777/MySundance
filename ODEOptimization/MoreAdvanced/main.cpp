#include "Sundance.hpp"
//#include "PlayaAmesosSolver.hpp"
#include "PlayaDenseSerialMatrix.hpp"
//#include "PlayaSerialVectorType.hpp"
#include "VientoSnapshotIO.hpp"
#include "denseSerialMatrixIO.hpp"


#include "integralOperator.hpp"
#include "ODECO.hpp"
#include "KKT_Transient_Channel.hpp"
#include "sensorData.hpp"

// Local Files
#include "MathematicaConverter.hpp"


// Standard Library Functions
using std::cout;
using std::endl;
using std::vector;

/**
 * getPoint(Point P) returns a CellFilter that will return the cell containing
 * the P(x,y) point
 */
CellFilter getPoint(Point P)
{
  CellFilter vertices = new DimensionalCellFilter(0);
  CellFilter xPos = vertices.coordSubset(0, P[0]);
  CellFilter yPos = vertices.coordSubset(1, P[1]);
  CellFilter vertex = xPos.intersection(yPos);

  return vertex;
}



    
int main(int argc, char *argv[]) 
{
  try
    {
      int verb = 0;      
      Sundance::setOption("verbosity", verb, "verbosity level");
      
      int nx = 25;
      Sundance::setOption("nx", nx, "Number of elements along each axis");

      int nSteps = 25;
      Sundance::setOption("nSteps", nSteps, "Number of time steps");

      //double gamma = 0.1;
      //Sundance::setOption("gamma", gamma, "perturbation to target");

      int quadOrder = 2;
      Sundance::setOption("quadOrder", quadOrder, "Order for the Gaussian Quadrature rule");

      // int Ru = 2;
      // Sundance::setOption("Ru", Ru, "Number of POD basis functions");

      double tFinal = 1.0;
      Sundance::setOption("tFinal",tFinal,"Final time value");

      int meshVerb = 0;
      Sundance::setOption("meshVerb",meshVerb,"Mesh verbosity level");

      // Probably will return eta to its previous location once we determine a good value
      double eta = 0.05;
      Sundance::setOption("eta",eta,"constant parameter on the regularization term");

      Sundance::init(&argc, &argv);

      // State the location of the necessary ROM files
      string ROM_base_dir = "/home/sirush/PhDResearch/MMS_Transient_Channel/Results/ROM/uRO/nx" + Teuchos::toString(nx) + "nt" + Teuchos::toString(nSteps);
      string POD_basis_fileprefix = ROM_base_dir + "/POD_basis";
      string uB_fileprefix = ROM_base_dir + "/uB";

      // Create the spatial mesh
      MeshType spatialMeshType = new Sundance::BasicSimplicialMeshType();
      double xmin = 0.0;
      double xmax = 1.0;
      double ymin = 0.0;
      double ymax = 1.0;
      MeshSource spatialMesher = new Sundance::PartitionedRectangleMesher(xmin,xmax,nx,1,ymin,ymax,nx,1,spatialMeshType,meshVerb);
      Mesh spatialMesh = spatialMesher.getMesh();

            // Define the time mesh
      MeshType timeMeshType = new BasicSimplicialMeshType();
      MeshSource timeMesher = new PartitionedLineMesher(0.0, tFinal, nSteps, timeMeshType);
      Mesh timeMesh = timeMesher.getMesh();

      // Set up appropriate cell filters
      CellFilter interior = new MaximalCellFilter();
      CellFilter verts = new DimensionalCellFilter(0);
      CellFilter left = verts.coordSubset(0,0.0);
      CellFilter right = verts.coordSubset(0,tFinal);

      VectorType<double> epetraVecType = new EpetraVectorType();

      
      // Read in b(t) from the ROM code for the MMS
      string b_filename = ROM_base_dir + "/b.txt";
      std::ifstream b_reader(b_filename, std::ios::in);
      TEUCHOS_TEST_FOR_EXCEPTION(!b_reader, std::runtime_error, "could not open file "
				 << b_filename << " for reading");

      int basisOrder;
      b_reader >> basisOrder;
      int numOfElements; // This is Ru
      b_reader >> numOfElements;
      int Ru = numOfElements;
      int numOfVectors; // equal to nSteps+1
      b_reader >> numOfVectors;
      string alias;
      b_reader >> alias;

      TEUCHOS_TEST_FOR_EXCEPTION(basisOrder != 1, std::runtime_error, "The order of the basis functions for time-dependent only functions is not Lagrange(1)");

      BasisFamily time_basis = new Lagrange(basisOrder);
      Array<BasisFamily> time_basisArray(Ru);
      for(int i = 0; i < Ru; i++)
	time_basisArray[i] = time_basis;
      DiscreteSpace time_DS(timeMesh, time_basisArray, epetraVecType);
      Expr b = new DiscreteFunction(time_DS, 0.0, alias);
      

      // BasisFamily b_bas = new Lagrange(basisOrder);
      // Array<BasisFamily> b_basArray(numOfElements);
      // for(int i = 0; i < numOfElements; i++)
      // 	b_basArray[i] = b_bas;
      // DiscreteSpace b_DS(timeMesh, b_basArray, epetraVecType);
      // Expr b = new DiscreteFunction(b_DS, 0.0, alias);

      // The get() here must somehow link to the DiscreteFunction's actual vector (shallow)
      Vector<double> b_Vec = getDiscreteFunctionVector(b);
      for(int i = 0; i < b_Vec.space().numLocalElements(); i++)
	b_reader >> b_Vec[i];

      // Read in the POD basis functions and uB from file
      Array<Expr> phi;
      phi.resize(Ru);
      for(int i = 0; i < Ru; i++)
	{
	  phi[i] = readSnap(POD_basis_fileprefix, i, spatialMesh);
	}
      Expr uB = readSnap(uB_fileprefix,0,spatialMesh);
      
     

      // Define t and its derivative
      Expr dt = new Derivative(0);
      //Expr t = new CoordExpr(0);

      /* Currently not using serial vectors anywhere
      Playa::VectorType<double> serialVecType = new Playa::SerialVectorType();
      VectorSpace<double> serialVecSpace = serialVecType.createEvenlyPartitionedSpace(MPIComm::self(),2);
      */
      
      string matrixAndTensorDir = "/home/sirush/PhDResearch/MMS_Transient_Channel/A_and_T/";
      string matrixFilename = matrixAndTensorDir + "A_nx" + Teuchos::toString(nx) + "nt" + Teuchos::toString(nSteps) + ".txt";

      ifstream is(matrixFilename, std::ios::in);
      TEUCHOS_TEST_FOR_EXCEPTION(!is, std::runtime_error, "could not open file " << matrixFilename << " for reading");

      // The information is written to file by denseSerialMatrixIO::writeDenseSerialMatrix
      // which does it row-wise
      Expr A;
      double value;
      for(int i = 0; i < Ru; i++)
	{
	  Expr row;
	  for(int j = 0; j < Ru; j++)
	    {
	      is >> value;
	      row.append(value);
	    }
	  A.append(row);
	}

      //cout << "Here is A: " << A << endl;

      // Form the transpose of A
      Expr At;
      for(int i = 0; i < Ru; i++)
	{
	  Expr row;
	  for(int j = 0; j < Ru; j++)
	    {
	      Expr tempValue = A[j][i];
	      row.append(tempValue);
	    }
	  At.append(row);
	}

      //cout << "Here is At: " << At << endl;

      string tensorPrefix = "T";
      Expr T;

      for(int k = 0; k < Ru; k++)
	{
	  // cout << "Starting the loop for k = " << k << endl;
	  string tensorFilename = matrixAndTensorDir + tensorPrefix + "[" + Teuchos::toString(k) + "]_nx" + Teuchos::toString(nx) + "nt" + Teuchos::toString(nSteps) + ".txt";
	  ifstream isTensor(tensorFilename, std::ios::in);
	  TEUCHOS_TEST_FOR_EXCEPTION(!isTensor, std::runtime_error, "could not open file " << tensorFilename << " for reading");

	  Expr Tk;
	  
	  for(int i = 0; i < Ru; i++)
	    {
	      Expr row;
	      for(int j = 0; j < Ru; j++)
		{
		  isTensor >> value;
		  row.append(value);
		}
	      Tk.append(row);
	    }
	  //cout << "Here is T[" << k << "]: " << Tk << endl;
	  T.append(Tk);
	}
      

      //cout << "Here is T: " << T << endl;




      
      // Read in alphaROM(t) from the ROM code for the MMS compartmentalize this in its own function
      // move the check from the end to here
      string alphaROM_filename = ROM_base_dir + "/alphaROM.txt";
      std::ifstream alphaROM_reader(alphaROM_filename, std::ios::in);
      TEUCHOS_TEST_FOR_EXCEPTION(!alphaROM_reader, std::runtime_error, "could not open file"
				 << alphaROM_filename << " for reading");

      alphaROM_reader >> basisOrder;
      alphaROM_reader >> numOfElements; // This is Ru
      alphaROM_reader >> numOfVectors;
      alphaROM_reader >> alias;


      TEUCHOS_TEST_FOR_EXCEPTION(basisOrder != 1, std::runtime_error, "The order of the basis functions for alphaROM  is not Lagrange(1)");
      TEUCHOS_TEST_FOR_EXCEPTION(numOfElements != Ru, std::runtime_error, "The value for Ru for alphaROM does not match that for b");

      Expr alphaROM = new DiscreteFunction(time_DS, 0.0, alias);

      // BasisFamily alphaROM_basis = new Lagrange(basisOrder);
      // Array<BasisFamily> alphaROM_basArray(numOfElements);
      // for(int i = 0; i < numOfElements; i++)
      // 	alphaROM_basArray[i] = alphaROM_basis;
      // DiscreteSpace alphaROM_DS(timeMesh, alphaROM_basArray, epetraVecType);
      // Expr alphaROM = new DiscreteFunction(alphaROM_DS, 0.0, alias);

      // The get() links to the DiscreteFunction's actual vector (shallow)
      Vector<double> alphaROM_Vec = getDiscreteFunctionVector(alphaROM);
      for(int i = 0; i < alphaROM_Vec.space().numLocalElements(); i++)
	alphaROM_reader >> alphaROM_Vec[i];

      



      // Create vstar (Measurement values)
      // Create an array that holds the * position values of the sensors
      Array<double> positionValues = tuple(0.2, 0.4, 0.6, 0.8);
      Array<Point> positionArray;
      // spatialDim() returns n for nD
      positionArray.resize(pow(positionValues.length(),spatialMesh.spatialDim() ));
      int Ns = positionArray.length(); // Number of sensors
      // Assigns values for (x1,y1) through (x1,yN), then goes to x2 and so forth
      for(int i=0; i<positionValues.length(); i++)
	for(int j=0; j<positionValues.length(); j++)
	  {
	    positionArray[j+positionValues.length()*i].resize(spatialMesh.spatialDim() );
	    positionArray[j+positionValues.length()*i][0] = positionValues[i];
	    positionArray[j+positionValues.length()*i][1] = positionValues[j];	  
	  }

      // Create the east unit vector
      double eastAngle = 0.0;
      Expr eastVec = List(cos(eastAngle), sin(eastAngle));

      // Read a snapshot for t_i into an Expr, then dot it with eastVec
      string snapshotDataDir = "/home/sirush/PhDResearch/MMS_Transient_Channel/Results/ForwardProblem/R=1forward_problem_TransientChannel_nx" + Teuchos::toString(nx) + "-nt-" + Teuchos::toString(nSteps) + "/";
      string filenamePrefix = "st-v";
      string snapshotFilenamePrefix = snapshotDataDir + filenamePrefix;


      // Starting here is in sensorData
      /* 
       * This check ensures that the sensor locations are a subset of the spatial mesh vertices
       */
      int numOfVertices = spatialMesh.numCells(0);
      for(int i = 0; i < Ns; i++)
	{
	  bool flag = false;
	  for(int j = 0; j < numOfVertices; j++)
	    {
	      if(positionArray[i][0] == spatialMesh.nodePosition(j)[0] &&
		 positionArray[i][1] == spatialMesh.nodePosition(j)[1] )
		{
		  flag = true;
		  break;
		}
	    }
	  TEUCHOS_TEST_FOR_EXCEPTION(!flag, std::runtime_error, "sensor locations are not a subset of spatial mesh vertices");	    
	}


      // this gets uForward(t_i=timeStep,x) for t_i = 0:tFinal
      Array<Expr> velocityDF;
      velocityDF.resize(nSteps+1);
      for(int timeStep = 0; timeStep < nSteps+1; timeStep++)
	{
	  velocityDF[timeStep] = readSnap(snapshotFilenamePrefix, timeStep, spatialMesh);
	}

      // Calculate the integral in the east direction to establish the Array of sensor values
      QuadratureFamily quad = new GaussianQuadrature(quadOrder);
      CellFilter pointFilter;
      Expr vstar;
      integralOperator S(spatialMesh, eastVec, quad);

      // Create the DiscreteSpace for vi
      //BasisFamily v_bas = new Lagrange(1);
      DiscreteSpace v_DS(timeMesh, time_basis, epetraVecType);

      
      // This takes the values from the forward simulation and stores them as the measurement
      // values at discrete time steps
      // The i loop is over the sensor locations, and the j loop is over the time steps
      for(int i = 0; i < Ns; i++) // Going through each vstar[i] to make it a DiscreteFunction
	{
	  Vector<double> values = v_DS.createVector();
	  pointFilter = getPoint(positionArray[i]);
	  
	  for(int j = 0; j < nSteps+1; j++) // Get the measurement values for vstar[i](t_j)
	    {
	      values[j] = S.staticDetect(pointFilter, velocityDF[j]);
      	      //cout << "The value of the integral eastVec*velocityDF over the point " << positionArray[i] << " at time step " << j << "  is " << values[j] << endl;	      
	    }
	  Expr vi = new DiscreteFunction(v_DS, values,"v["+Teuchos::toString(i)+"]");
	  vstar.append(vi);
	}

      // Here ends the sensorData file



      // How to replicate this?
      // /* "exact solution" produce by Mathematica's FindMinimum. Accurate to
      //  * about 1.0e-6. */
      // Array<double> alphaExact = Teuchos::tuple(
      // 						0.5000000289639744 + gamma*
      // 						(0.064717637614702 + 
      // 						 (-0.00044913462946438883 - 
      // 						  0.00001613396940867341*gamma)*gamma),
      // 						0.3333333245793077 + gamma*
      // 						(-0.10739738478174649 + 
      // 						 (-0.0015852504354095308 - 
      // 						  0.000019917356023471516*gamma)*gamma)
      // 						);




      // Form the KKT system
      // State Variable: alpha
      // Adjoint Variable: lambda
      // Design Variable: p

      // From here is in in KKTBase
      Expr alpha, lambda, p;
      for(int i=0; i<Ru; i++)
      	{
      	  alpha.append(new UnknownFunction(time_basis, "alpha_"+Teuchos::toString(i)));
      	  lambda.append(new UnknownFunction(time_basis, "lambda"));
      	  p.append(new UnknownFunction(time_basis, "p"));
      	}
      
      Expr alphaHat, lambdaHat, pHat;
      for(int i=0; i<Ru; i++)
      	{
      	  alphaHat.append(new TestFunction(time_basis));
      	  lambdaHat.append(new TestFunction(time_basis));
      	  pHat.append(new TestFunction(time_basis));
      	}

      Array<BasisFamily> ODECO_basisArray(3*Ru);
      for (int i=0; i<ODECO_basisArray.size(); i++)
	ODECO_basisArray[i]=time_basis;
      DiscreteSpace ODECO_DS(timeMesh, ODECO_basisArray, epetraVecType);
      // To here is in KKTBase


      
      Expr U0 = new DiscreteFunction(ODECO_DS, 0.0);


      /*
       * The derivation of the state, adjoint, and design equations must be done by hand on a per-problem basis :(
       */

      /* state equation & BCs */
      // H(T,alpha) is implemented in sundance as T*alpha*alpha
      Expr stateEqn = Integral(interior, lambdaHat*(dt*alpha - A*alpha - T*alpha*alpha -b), quad);
      Expr stateBC = EssentialBC(left, lambdaHat*(alpha-p), quad);

      /* adjoint equation & BCs, derived by hand */
      //Expr adjointEqn = Integral(interior, alphaHat*(alpha-vstar) - alphaHat*(dt*lambda+At*lambda) - lambda*(T*alphaHat*alpha + T*alpha*alphaHat), quad);
      //Expr adjointBC = EssentialBC(right, alphaHat*lambda, quad);

      /* Store the POD basis functions as columns of a "matrix"
       * Currently I am not using this anywhere, but I am keeping it as an example
       * of how to write the POD basis functions as columns of a matrix
      Expr PHI;
      for(int dim = 0; dim < 2; dim++)
	{
	  Expr row;
	  for(int r = 0; r < Ru; r++)
	    {
	      row.append(phi[r][dim]);
	    }
	  PHI.append(row);
	}
      */
      
      // S_i(u) = S_i(uB) + sum( alpha[r] S_i( phi[r] ) )
      Expr Su;
      //      integralOperator S(spatialMesh, eastVec, quad);
      for(int i = 0; i < Ns; i++)
	{
	  // Calculate S_i(uB)
	  CellFilter location = getPoint(positionArray[i]);
	  Expr Si = S.staticDetect(location, uB);
	  for(int r = 0; r < Ru; r++)
	    {
	      Si = Si + alpha[r]*S.staticDetect(location,phi[r]);
	    }
	  Su.append(Si);
	}

      // The purpose of this section of code is to calculate (Su - vstar, S(PHI*alphaHat) )
      // For the moment let the IP be the dot product
      Expr SPHIaHat;
      for(int i = 0; i < Ns; i++)
	{
	  CellFilter location = getPoint(positionArray[i]);
	  Expr Si = 0.0;
	  for(int r = 0; r < Ru; r++)
	    {
	      Si = Si + alphaHat[r]*S.staticDetect(location,phi[r]);
	    }
	  SPHIaHat.append(Si);
	}
      
      //      Expr adjointEqn = Integral(interior, (Su - vstar)*SPHIaHat - ( alphaHat*(dt*lambda + At*lambda) + lambda*(T*alphaHat*alpha + T*alpha*alphaHat) ), quad);
      //Might have made an error on the signs
      Expr adjointEqn = Integral(interior, (Su - vstar)*SPHIaHat - ( alphaHat*(dt*lambda - At*lambda) - lambda*(T*alphaHat*alpha + T*alpha*alphaHat) ), quad);
      //double eta = 0.05; // Need to ask about this value
      Expr regTerm = 0.0;
      // summation over i
      for(int i = 0; i < Ru; i++)
	{
	  FunctionalEvaluator temp(spatialMesh, Integral(interior, uB*phi[i], quad));
	  regTerm = regTerm + alphaHat[i]*temp.evaluate();
	}

      // // summation over j and k
      // for(int j = 0; j < Ru; j++)
      // 	{
      // 	  for(int k = 0; k < Ru; k++)
      // 	    {
      // 	      FunctionalEvaluator temp(spatialMesh, Integral(interior, phi[j]*phi[k],quad));
      // 	      cout << "(phi[" << j << "], phi[" << k << "]) = " << temp.evaluate() << endl;
      // 	      regTerm = regTerm + alpha[j]*alphaHat[k]*temp.evaluate();
      // 	    }
      // 	}

      // Summation over j and k simplifying due to orthogonality of the basis functions
      for(int j = 0; j < Ru; j++)
	{
	  regTerm = regTerm + alpha[j]*alphaHat[j];
	}

      // Multiply by the constant term eta
      regTerm = regTerm*eta;
      adjointEqn = adjointEqn + Integral(interior, regTerm, quad);
      Expr adjointBC = EssentialBC(right, alphaHat*lambda, quad);



      
      /* design equation and BC */
      /* -- the (alpha')*(alphaHat') term enforces constancy of alpha in time */
      // currently adding only 0.5 int( dp/dt * dp/dt ) without eta
      Expr designEqn = Integral(interior, (dt*p)*(dt*pHat), quad);
      Expr designBC = EssentialBC(left, pHat*lambda, quad);

      /* Combine the equations and BCs to form the KKT system  */
      Expr eqn = stateEqn + adjointEqn + designEqn;
      Expr bc = stateBC + adjointBC + designBC;

      /* create the NLP */
      NonlinearProblem NLP(timeMesh, eqn, bc, List(lambdaHat, alphaHat, pHat),
      			   List(alpha, lambda, p), U0, epetraVecType);

      NonlinearSolver<double> solver 
        = NonlinearSolverBuilder::createSolver("playa-newton-amesos.xml");

      SolverState<double> state = NLP.solve(solver);
      

      TEUCHOS_TEST_FOR_EXCEPTION(state.finalState() != SolveConverged,
                                 std::runtime_error,
                                 "Nonlinear solve failed to converge: message="
                                 << state.finalMsg());

      // The solve for alpha values will be the first Ru components of U0
      // We need to now build the approximation to the velocity uOpt
      Expr uOpt = uB;
      for(int r = 0; r < Ru; r++)
	{
	  uOpt = uOpt + U0[r]*phi[r];
	}

      /* 
       * For the moment, don't know how to compare to uEx and pEx
       *
       // Compare to uExact ask how to do this
      Expr x = new CoordExpr(0,"x");
      Expr y = new CoordExpr(1,"y");
      Expr uExact = List(1 - (Pi*(-1 + Power(x,2))*(8*Sin((Pi*t)/2.)*Power(Sin(Pi*x),2)*Sin(4*Pi*y) +45*Sin(Pi*t)*Power(Sin(2*Pi*x),2)*Sin(6*Pi*y)))/200.,(4*Sin((Pi*t)/2.)*Sin(Pi*x)*(Pi*(-1 + Power(x,2))*Cos(Pi*x) + x*Sin(Pi*x))*Power(Sin(2*Pi*y),2) +15*Sin(Pi*t)*Sin(2*Pi*x)*(2*Pi*(-1 + Power(x,2))*Cos(2*Pi*x) + x*Sin(2*Pi*x))*Power(Sin(3*Pi*y),2))/100.); 

      Expr pExact = 0.0;
      */
      
      // cout << "t-x = " << t-x << endl; // so yeah, sundance thinks t and x are the same

      /*
       * Since uOpt and uRO only differe in their values for alpha, we will compare those
       * alphaOPT = U0
       */
      Expr alphaOPT = U0[0];
      for(int r=1; r < Ru; r++)
	alphaOPT.append(U0[r]);

      double alphaError = L2Norm(timeMesh, interior, alphaROM-alphaOPT, quad);
      cout << "||alphaROM - alphaOPT||_2 = " << alphaError << endl;

      cout << "Just checking against an expression with tiny values (0.001, 0.0001) " << endl;
      Expr tester = List(0.001, 0.0001);
      cout << "||alphaROM - tester||_2 = " << L2Norm(timeMesh, interior, alphaROM - tester, quad) << endl;

      
      // Array<double> alphaNum = Teuchos::tuple(
      // 					      L2Norm(mesh, left, U0[0], quad),
      // 					      L2Norm(mesh, left, U0[1], quad)
      // 					      );
      // for (int j=0; j<2; j++)
      // 	{
      // 	  Tabs tab1;
      // 	  Out::os() << tab1 << "Alpha[" << j << "]: exact="
      // 		    << alphaExact[j]
      // 		    << ", numerical " << alphaNum[j]
      // 		    << ", error=" << fabs(alphaExact[j] - alphaNum[j])
      // 		    << endl;
      // 	}
      
      
      // FieldWriter writer = new DSVWriter("2DLinearOpt-.dat");
      // writer.addMesh(mesh);
      // writer.addField("x[0]", new ExprFieldWrapper(U0[0]));
      // writer.addField("x[1]", new ExprFieldWrapper(U0[1]));
      // writer.addField("lambda[0]", new ExprFieldWrapper(U0[2]));
      // writer.addField("lambda[1]", new ExprFieldWrapper(U0[3]));
      // writer.addField("alpha[0]", new ExprFieldWrapper(U0[4]));
      // writer.addField("alpha[1]", new ExprFieldWrapper(U0[5]));
      // writer.write();
      
      // Array<double> alphaNum = Teuchos::tuple(
      // 					      L2Norm(mesh, left, U0[0], quad),
      // 					      L2Norm(mesh, left, U0[1], quad)
      // 					      );
      // for (int j=0; j<2; j++)
      // 	{
      // 	  Tabs tab1;
      // 	  Out::os() << tab1 << "Alpha[" << j << "]: exact="
      // 		    << alphaExact[j]
      // 		    << ", numerical " << alphaNum[j]
      // 		    << ", error=" << fabs(alphaExact[j] - alphaNum[j])
      // 		    << endl;
      // 	}
    }
  catch(std::exception& ex)
    {
      cerr << "main() caught exception: " << ex.what() << endl;
    }
  Sundance::finalize();
}

    

    
    





