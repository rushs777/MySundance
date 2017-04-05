#include "MyNLO.hpp"

MyNLO::MyNLO(MMSQuadODE f, double h) : NonlinearOperatorBase(f.space(), f.space()), f_(f) 
  {
    uPrev_ = f.space().createMember();
    setEvalPt(getInitialGuess());

    h_ = h;
    tPrev_ = 0.0;
    tNext_ = h;
  }

void MyNLO::set_tPrev(const double t)
  {
    tPrev_ = t;
    tNext_ = tPrev_ + h_;
  }

LinearOperator<double> MyNLO::computeJacobianAndFunction(Vector<double>& functionValue) const
  {
    LinearOperator<double> Jf;
    // Calculate f(t_n, u_n)
    Vector<double> fPrev = f_.eval1(tPrev_, uPrev_, Jf);
    // Calculate f(t_{n+1}, x^n)
    Vector<double> fNext = f_.eval1(tNext_, currentEvalPt(), Jf); 
    // Calculate F(x^n)
    functionValue = currentEvalPt() - uPrev_ - (h_/2.0)*(fPrev + fNext);
    // Calculate JF(x^n) = I - (h/2)Jf(x^n)
    LinearOperator<double> JF(rcp(new DenseSerialMatrix(f_.space(), f_.space())));
    RCP<DenseSerialMatrix> JFptr = DenseSerialMatrix::getConcretePtr(JF);
    RCP<DenseSerialMatrix> Jfptr = DenseSerialMatrix::getConcretePtr(Jf);

    for(int i = 0; i<Jfptr->numRows(); i++)
      for(int j = 0; j<Jfptr->numCols(); j++)
	{
	  if(i==j)
	    JFptr->setElement(i,i, 1.0 - (h_/2.0)*Jfptr->getElement(i,i));
	  if(i!=j)
	    JFptr->setElement(i,j, (-h_/2.0)*Jfptr->getElement(i,j));
	}

    return JF;
  }
