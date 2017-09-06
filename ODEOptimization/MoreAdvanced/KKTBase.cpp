#include "KKTBase.hpp"


KKTBase::KKTBase(string ROM_base_dir, Mesh spatialMesh, Mesh timeMesh, double tFinal, int verbosity)
  :  ROM_base_dir_(ROM_base_dir),
     spatialMesh_(spatialMesh),
     timeMesh_(timeMesh),
     tFinal_(tFinal),
     verbosity_(verbosity)
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
}

void KKTBase::initialize_b()
{
  string b_filename = ROM_base_dir_ + "/b.txt";

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
  string POD_basis_fileprefix = ROM_base_dir_ + "/POD_basis";

  phi_.resize(Ru_);
  for(int i = 0; i < Ru_; i++)
    {
      phi_[i] = readSnap(POD_basis_fileprefix, i, spatialMesh_);
    }
}


void KKTBase::initialize_uB()
{
  string uB_fileprefix = ROM_base_dir_ + "/uB";
  uB_ = readSnap(uB_fileprefix,0,spatialMesh_);
}


void KKTBase::initialize_A_and_T()
{
  string matrixFilename = ROM_base_dir_ + "/A.txt";

  // ifstream is(matrixFilename, std::ios::in);
  // TEUCHOS_TEST_FOR_EXCEPTION(!is, std::runtime_error, "could not open file " << matrixFilename << " for reading");

  // // The information is written to file by denseSerialMatrixIO::writeDenseSerialMatrix
  // // which does it row-wise
  // // Add this code to denseSerialMatrixIO.cpp
  // double value;
  // for(int i = 0; i < Ru_; i++)
  //   {
  //     Expr row;
  //     for(int j = 0; j < Ru_; j++)
  // 	{
  // 	  is >> value;
  // 	  row.append(value);
  // 	}
  //     A_.append(row);
  //   }

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

  //double value;
  for(int k = 0; k < Ru_; k++)
    {
      // cout << "Starting the loop for k = " << k << endl;
      string tensorFilename = ROM_base_dir_ + "/T[" + Teuchos::toString(k) + "].txt";
      // ifstream isTensor(tensorFilename, std::ios::in);
      // TEUCHOS_TEST_FOR_EXCEPTION(!isTensor, std::runtime_error, "could not open file " << tensorFilename << " for reading");

      Expr Tk;
	  
      // for(int i = 0; i < Ru_; i++)
      // 	{
      // 	  Expr row;
      // 	  for(int j = 0; j < Ru_; j++)
      // 	    {
      // 	      isTensor >> value;
      // 	      row.append(value);
      // 	    }
      // 	  Tk.append(row);
      // 	}

      readListExprMatrix(Tk, Ru_, tensorFilename);
      
      SUNDANCE_ROOT_MSG2(verbosity_, "T[" << k << "]: " << endl << Tk << endl);
      T_.append(Tk);
    }
      
  SUNDANCE_ROOT_MSG2(verbosity_, "T: " << endl << T_ << endl);
}

