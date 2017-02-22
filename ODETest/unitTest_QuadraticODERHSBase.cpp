#include "Teuchos_GlobalMPISession.hpp"
#include "QuadraticODERHSBase.hpp"

#include <vector>
using std::vector;

/*
using namespace Teuchos;
using namespace Playa;
using namespace PlayaExprTemplates;
*/

class TestQuadODE : public QuadraticODERHSBase
{
public:
  TestQuadODE() : QuadraticODERHSBase(3, 2) {}

  Vector<double> evalForceTerm(const double& t) const
  {
    Vector<double> rtn = space().createMember();
    rtn[0] = 1 + t*(-6 + t*(3 + (-1 - 3*t)*t));
    rtn[1] = 4 + t*(4 + t*(-8 + t*(-2 + 2*t)));
    rtn[2] = -2 + t*(-3 + t*(-1 + (-3 - 2*t)*t));
    
    return rtn*0.0;
  }

  void fillMatrixAndTensor(RCP<DenseSerialMatrix>& A,
			   Array<RCP<DenseSerialMatrix> >& T)
  {
    vector<vector<double> > AVals {{2,-1,0},{-1,2,-1},{0,-1,2}};
    vector<vector<vector<double> > > TVals = {{{-3,-4,4},{-3,1,0},{3,1,3}},
					      {{-2,3,4},{2,-2,3},{2,2,4}},
					      {{2,-3,4},{-1,-1,3},{-2,-1,-1}}};
    
    A->fill(AVals);
    for (int i=0; i<AVals.size(); i++)
      {
	T[i]->fill(TVals[i]);
      }
  }
};



int main(int argc, char *argv[]) 
{
  try
  {
    GlobalMPISession session(&argc, &argv);

    TestQuadODE f;
    f.initialize();

    double t = 1.0;
    
    Vector<double> u = f.space().createMember();
    u[0] = 1.0;
    u[1] = 2.0;
    u[2] = 3.0;

    LinearOperator<double> J;
    Vector<double> fVal = f.eval1(t, u, J);

    Out::root() << "f = " << endl << fVal << endl;
    Out::root() << "J = " << endl << J << endl;
  }
  catch(std::exception& e)
  {
    Out::root() << "Caught exception: " << e.what() << std::endl;
    return -1;
  }
}

