#include "KKTBase.hpp"

KKTBase::KKTBase(Mesh timeMesh, int Ru) :
  Ru_(Ru),
  interior_(new MaximalCellFilter())
{
  BasisFamily time_basis = new Lagrange(1);
  VectorType<double> epetraVecType = new EpetraVectorType();
  
  for(int i=0; i<Ru_; i++)
    {
      alpha_.append(new UnknownFunction(time_basis, "alpha_"+Teuchos::toString(i)));
      lambda_.append(new UnknownFunction(time_basis, "lambda"));
      p_.append(new UnknownFunction(time_basis, "p"));
    }
  
  for(int i=0; i<Ru; i++)
    {
      alphaHat_.append(new TestFunction(time_basis));
      lambdaHat_.append(new TestFunction(time_basis));
      pHat_.append(new TestFunction(time_basis));
    }

  Array<BasisFamily> ODECO_basisArray(3*Ru_);
  for (int i=0; i<ODECO_basisArray.size(); i++)
    ODECO_basisArray[i]=time_basis;
  DiscreteSpace ODECO_DS_(timeMesh, ODECO_basisArray, epetraVecType);
  //Expr U0 = new DiscreteFunction(ODECO_DS, 0.0);

  Expr dt_ = new Derivative(0);
}
