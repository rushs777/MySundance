#include "KKTBase.hpp"


KKTBase::KKTBase(string POD_DataDir, Mesh spatialMesh, Mesh timeMesh, double tFinal, int nSteps, int verbosity)
  :  POD_DataDir_(POD_DataDir),
     spatialMesh_(spatialMesh),
     timeMesh_(timeMesh),
     tFinal_(tFinal),
     nSteps_(nSteps),
     verbosity_(verbosity),
     uOPT_()
{
  dt_ = new Derivative(0);
  interior_ = new MaximalCellFilter();
}


void KKTBase::initialize_vars()
{ 
  for(int i=0; i<Ru_; i++)
    {
      alpha_.append(new UnknownFunction(time_basis_, "alpha_"+Teuchos::toString(i)));
      lambda_.append(new UnknownFunction(time_basis_, "lambda"));
      p_.append(new UnknownFunction(time_basis_, "p"));
    }
  
  for(int i=0; i<Ru_; i++)
    {
      alphaHat_.append(new TestFunction(time_basis_));
      lambdaHat_.append(new TestFunction(time_basis_));
      pHat_.append(new TestFunction(time_basis_));
    }
  
}


void KKTBase::initialize()
{
  epetraVecType_ = new EpetraVectorType();
  initialize_b();
  initialize_phi();
  initialize_uB();
  initialize_A_and_T();
  initialize_vars();
  
  for(int i=0; i<nSteps_+1; i++)
    uOPT_.append(uB_);
}

void KKTBase::initialize_b()
{
  string b_filename = POD_DataDir_ + "b.txt";

  // Read in the preamble information for b(t) from the ROM code
  std::ifstream b_reader(b_filename, std::ios::in);
  TEUCHOS_TEST_FOR_EXCEPTION(!b_reader, std::runtime_error, "could not open file "
			     << b_filename << " for reading");
  int basisOrder;
  b_reader >> basisOrder;
  b_reader >> Ru_; // Ru is the number of elements in vector-valued function b
  int numOfVectors; // equal to nSteps+1
  b_reader >> numOfVectors;
  string alias;
  b_reader >> alias;

  // Check to make sure that the time-dependent only functions are using Lagrange(1)
  TEUCHOS_TEST_FOR_EXCEPTION(basisOrder != 1, std::runtime_error, "The order of the basis functions for time-dependent only functions is not Lagrange(1)");


  // Create the time DiscreteSpace
  time_basis_ = new Lagrange(basisOrder);
  Array<BasisFamily> time_basisArray(Ru_);
  for(int i = 0; i < Ru_; i++)
    time_basisArray[i] = time_basis_;
  //  time_DS_ = new DiscreteSpace(timeMesh_, time_basisArray, epetraVecType_);
  DiscreteSpace time_DS(timeMesh_, time_basisArray, epetraVecType_);
  
  // Establish the values for b(t)
  b_ = new DiscreteFunction(time_DS, 0.0, alias);
  Vector<double> b_Vec = getDiscreteFunctionVector(b_);
  for(int i = 0; i < b_Vec.space().numLocalElements(); i++)
    b_reader >> b_Vec[i];  

}


void KKTBase::initialize_phi()
{
  // Read in the POD from file
  SUNDANCE_ROOT_MSG1(verbosity_, "Reading in the reduced-order basis for velocity");
  string POD_basis_fileprefix = POD_DataDir_ + "POD_basis";

  Ru_ = 10000;
  for(int r = 0; r < Ru_; r++)
    {
      try
	{
	  phi_.push_back( readSnap(POD_basis_fileprefix, r, spatialMesh_ ) );
	}
      catch (std::runtime_error& e)
	{
	  Ru_ = r;
	}
    }
  
  SUNDANCE_ROOT_MSG1(verbosity_, "Found " << Ru_ << " reduced-order basis functions for the given resoultion and tolerance"); 
}


void KKTBase::initialize_uB()
{
  string uB_fileprefix = POD_DataDir_ + "uB";
  uB_ = readSnap(uB_fileprefix,0,spatialMesh_);
}


void KKTBase::initialize_A_and_T()
{
  string matrixFilename = POD_DataDir_ + "A.txt";

  readListExprMatrix(A_, Ru_, matrixFilename);

  SUNDANCE_ROOT_MSG2(verbosity_, "A: " << endl << A_ << endl);

  // Form the transpose of A
  for(int i = 0; i < Ru_; i++)
    {
      Expr row;
      for(int j = 0; j < Ru_; j++)
	{
	  Expr tempValue = A_[j][i];
	  row.append(tempValue);
	}
      At_.append(row);
    }

  SUNDANCE_ROOT_MSG2(verbosity_, "At: " << endl << At_ << endl);

  for(int k = 0; k < Ru_; k++)
    {
      string tensorFilename = POD_DataDir_ + "T[" + Teuchos::toString(k) + "].txt";
      Expr Tk;
      readListExprMatrix(Tk, Ru_, tensorFilename);
      SUNDANCE_ROOT_MSG2(verbosity_, "T[" << k << "]: " << endl << Tk << endl);
      T_.append(Tk);
    }
      
  SUNDANCE_ROOT_MSG2(verbosity_, "T: " << endl << T_ << endl);
}

