#include "QuadraticODERHSBase.hpp"

QuadraticODERHSBase::QuadraticODERHSBase(int n, int verbosity) 
: ODERHSBase(n), 
  A_(rcp(new DenseSerialMatrix(space(), space()))),
  T_(n), 
  initialized_(false),
  verbosity_(verbosity)
{}


void QuadraticODERHSBase::initialize()
{
    int n = space().dim();
    RCP<DenseSerialMatrix> APtr = DenseSerialMatrix::getConcretePtr(A_);
    Array<RCP<DenseSerialMatrix> > TPtr(n);
    for (int i=0; i<n; i++)
      {
	T_[i] = LinearOperator<double>(rcp(new DenseSerialMatrix(space(), space())));
	TPtr[i] = DenseSerialMatrix::getConcretePtr(T_[i]);
      }


    fillMatrixAndTensor(APtr, TPtr);

    if(verbosity_ >= 2)
    {
        Out::root() << "A=" << endl << A_ << endl;
        Out::root() << "APtr= " << endl;
        APtr->print(Out::root());
        for (int i=0; i<n; i++) Out::root() << "T[" << i << "] = " << endl << T_[i] << endl;
    }

    initialized_ = true;
}


Vector<double> QuadraticODERHSBase::eval1(const double& t,
		       const Vector<double>& u,
		       LinearOperator<double>& J) const
{
    TEUCHOS_TEST_FOR_EXCEPT(!initialized_);
    int n = space().dim();
    Vector<double> rtn =  A_*u + evalForceTerm(t);
    RCP<const DenseSerialMatrix> APtr = DenseSerialMatrix::getConcretePtr(A_);
    RCP<DenseSerialMatrix> JPtr = rcp(new DenseSerialMatrix(*APtr)); // deep copy of A
    J = LinearOperator<double>(JPtr);

    Array<Vector<double> > JCol(n);
    for (int i=0; i<n; i++)
      {
	rtn[i] += u*(T_[i]*u);
	JCol[i] = T_[i]*u + T_[i].transpose()*u;
      }

    for (int i=0; i<n; i++)
      {
	for (int j=0; j<n; j++)
	  {
	    double Jij = JPtr->getElement(i,j);
	    JPtr->setElement(i, j, Jij + JCol[j][i]); 
	  }
      }

    return rtn;
}

int QuadraticODERHSBase::getVerbosity() const 
{return verbosity_;}
