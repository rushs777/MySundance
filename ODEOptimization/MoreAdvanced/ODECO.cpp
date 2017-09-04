#include "ODECO.hpp"

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

ODECO::ODECO(string ROM_base_dir, string matrixAndTensorDir, Mesh spatialMesh, int verbosity)
  :  ROM_base_dir_(ROM_base_dir),
     spatialMesh_(spatialMesh),
     timeMesh_(timeMesh),
     matrixAndTensorDir_(matrixAndTensorDir),
     verbosity_(verbosity)
{


  // These values are depending on hard-coding in uRO.cpp


}

void ODECO::initialize()
{
  epetraVecType_ = new EpetraVectorType();
  initialize_b();
  initialize_phi();
  initialize_uB();
  initialize_A_and_T();
}

void ODECO::initialize_b()
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
  time_DS_ = new DiscreteSpace(timeMesh_, time_basisArray, epetraVecType_);

  // Establish the values for b(t)
  b_ = new DiscreteFunction(time_DS_, 0.0, alias);
  Vector<double> b_Vec = getDiscreteFunctionVector(b_);
  for(int i = 0; i < b_Vec.space().numLocalElements(); i++)
    b_reader >> b_Vec[i];  

}

void ODECO::initialize_phi()
{
  string POD_basis_fileprefix = ROM_base_dir_ + "/POD_basis";

  phi_.resize(Ru_);
  for(int i = 0; i < Ru_; i++)
    {
      phi_[i] = readSnap(POD_basis_fileprefix, i, spatialMesh_);
    }
}

void ODECO::initialize_uB()
{
  string uB_fileprefix = ROM_base_dir_ + "/uB";
  uB_ = readSnap(uB_fileprefix,0,spatialMesh_);
}

void ODECO::initialize_A_and_T()
{
  string matrixFilename = matrixAndTensorDir_ + "A_nx" + Teuchos::toString(nx) + "nt" + Teuchos::toString(nSteps) + ".txt";

  ifstream is(matrixFilename, std::ios::in);
  TEUCHOS_TEST_FOR_EXCEPTION(!is, std::runtime_error, "could not open file " << matrixFilename << " for reading");

  // The information is written to file by denseSerialMatrixIO::writeDenseSerialMatrix
  // which does it row-wise
  // Add this code to denseSerialMatrixIO.cpp
  double value;
  for(int i = 0; i < Ru_; i++)
    {
      Expr row;
      for(int j = 0; j < Ru_; j++)
	{
	  is >> value;
	  row.append(value);
	}
      A_.append(row);
    }

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
      // cout << "Starting the loop for k = " << k << endl;
      string tensorFilename = matrixAndTensorDir_ + "T[" + Teuchos::toString(k) + "]_nx" + Teuchos::toString(nx) + "nt" + Teuchos::toString(nSteps) + ".txt";
      ifstream isTensor(tensorFilename, std::ios::in);
      TEUCHOS_TEST_FOR_EXCEPTION(!isTensor, std::runtime_error, "could not open file " << tensorFilename << " for reading");

      Expr Tk;
	  
      for(int i = 0; i < Ru_; i++)
	{
	  Expr row;
	  for(int j = 0; j < Ru_; j++)
	    {
	      isTensor >> value;
	      row.append(value);
	    }
	  Tk.append(row);
	}
      SUNDANCE_ROOT_MSG2(verbosity_, "T[" << k << "]: " << endl << Tk << endl);
      T_.append(Tk);
    }
      
  SUNDANCE_ROOT_MSG2(verbosity_, "T: " << endl << T << endl);
}

Expr ODECO::get_b()
{
  return b_;
}

Array<Expr> ODECO::get_phi()
{
  return phi_;
}

Expr ODECO::get_uB()
{
  return uB_;
}

Expr ODECO::get_A()
{
  return A_;
}

Expr ODECO::get_At()
{
  return At_;
}

Expr ODECO::get_T()
{
  return T_;
}


double errorCheck(Expr alphaOPT)
{
  // Read in alphaROM(t) from the ROM code for the MMS compartmentalize this in its own function
  string alphaROM_filename = ROM_base_dir_ + "/alphaROM.txt";
  std::ifstream alphaROM_reader(alphaROM_filename, std::ios::in);
  TEUCHOS_TEST_FOR_EXCEPTION(!alphaROM_reader, std::runtime_error, "could not open file"
			     << alphaROM_filename << " for reading");

  int basisOrder;
  alphaROM_reader >> basisOrder;
  int numOfElements;
  alphaROM_reader >> numOfElements; // This is Ru
  int numOfVectors;
  alphaROM_reader >> numOfVectors;
  string alias;
  alphaROM_reader >> alias;


  TEUCHOS_TEST_FOR_EXCEPTION(basisOrder != 1, std::runtime_error, "The order of the basis functions for alphaROM  is not Lagrange(1)");
  TEUCHOS_TEST_FOR_EXCEPTION(numOfElements != Ru, std::runtime_error, "The value for Ru for alphaROM does not match that for b");

  Expr alphaROM = new DiscreteFunction(time_DS_, 0.0, alias);
  Vector<double> alphaROM_Vec = getDiscreteFunctionVector(alphaROM);
  for(int i = 0; i < alphaROM_Vec.space().numLocalElements(); i++)
    alphaROM_reader >> alphaROM_Vec[i];
}
