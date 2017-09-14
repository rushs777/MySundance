#include "Sundance.hpp"
#include "PlayaSerialVectorType.hpp"
#include "PlayaDenseSerialMatrix.hpp"
#include "PlayaDenseMatrixMatrixProduct.hpp"
//#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaSimpleAddedOpImpl.hpp"
// Add the Declaration; but not for adding

using std::cout;
using std::endl;
using std::vector;

int main(int argc, char *argv[])
{
  cout << "Hello, world" << endl;
  /*
   * Things I need: snapshot matrix W, DiscreteSpace ds for basis and mesh
   */

  Playa::VectorType<double> epetraVecType = new Playa::EpetraVectorType();
  Playa::VectorType<double> serialVecType = new Playa::SerialVectorType();

  
  /*
   * This section of the code creates a matrix A = W^T S W where S is
   * symmetric, positive definite
   * From it we got the SVD, A = U Sigma Vt, where U == V
   * and sigma >= 0
   */
  // Create the matrix A
  // Note: A = W^T S W, where the mass matrix S is symmetric, positive definite
  VectorSpace<double> cols = serialVecType.createEvenlyPartitionedSpace(MPIComm::world(), 3);
  VectorSpace<double> rows = serialVecType.createEvenlyPartitionedSpace(MPIComm::world(), 3);
  RCP<MatrixFactory<double> > mf = serialVecType.createMatrixFactory(cols, rows);
  LinearOperator<double> A = mf->createMatrix();
  RCP<DenseSerialMatrix> Aptr = DenseSerialMatrix::getConcretePtr(A);
  vector<vector<double> > AVals = {{26,38,50},{38,56,74},{50,74,98}};
  Aptr->fill(AVals);

  // Get the SVD of A
  LinearOperator<double> U;
  Vector<double> sigma;
  LinearOperator<double> Vt;
  denseSVD(A,U,sigma,Vt);

  cout << "U is : " << endl << U << endl;
  cout << "sigma is : " << endl << sigma << endl;
  cout << "Vt is : " << endl << Vt << endl;

  /*
   * Now we will create a LinearProblem to ensure that this matrix is SPD
   * Dr Long start here
   */
  // Define our mesh
  int nx = 1;
  MeshType meshType = new Sundance::BasicSimplicialMeshType();
  double xmin = 0.0;
  double xmax = 1.0;
  double ymin = 0.0;
  double ymax = 1.0;
  MeshSource mesher = new Sundance::PartitionedRectangleMesher(xmin,xmax,nx,1,
							       ymin,ymax,nx,1,
							       meshType,0);
  Mesh mesh = mesher.getMesh();

  // Create a BasisFamily to express our solution
  Array<Sundance::BasisFamily> basis = List(new Sundance::Lagrange(1),
					     new Sundance::Lagrange(1)); 
  // Define our vector type
  Sundance::DiscreteSpace ds(mesh, basis, epetraVecType);

  CellFilter interior = new MaximalCellFilter();

  // Test Functions
  Teuchos::Array<Expr> v;
  for(int i = 0; i < basis.size(); i++)
    {
      v.push_back(new TestFunction(basis[i], "v[" + Teuchos::toString(i) + "]"));
    }
	
  Sundance::Expr vlist = new Sundance::ListExpr(v);

  // Unknown Functions
  Teuchos::Array<Expr> u;
  for(int i = 0; i < basis.size(); i++)
    {
      u.push_back(new UnknownFunction(basis[i], "u[" + Teuchos::toString(i) + "]"));
    }

  Sundance::Expr ulist = new Sundance::ListExpr(u);

  // Evaluation scheme. The parameter says for what degree of polynomials it will be exact for
  Sundance::QuadratureFamily quad = new Sundance::GaussianQuadrature(4);
  // Define what you want to integrate
  Sundance::Expr integrand = vlist*ulist;
  Sundance::Expr eqn = Integral(interior,integrand,quad);
  // Define Empty BC since I want the mass matrix
  Sundance::Expr bc;
  // Define the problem
  Sundance::LinearProblem prob(mesh,eqn,bc,vlist,ulist,epetraVecType);

  // Get the mass matrix
  Playa::LinearOperator<double> S = prob.getOperator();

  cout << "S is " << S.range().dim() << " by " << S.domain().dim() << endl;
  cout << "S = " << S << endl;

  VectorSpace<double> Wcols = serialVecType.createEvenlyPartitionedSpace(MPIComm::world(), S.domain().dim() );
  VectorSpace<double> Wrows = serialVecType.createEvenlyPartitionedSpace(MPIComm::world(), S.range().dim() );
  RCP<MatrixFactory<double> > Wmf = serialVecType.createMatrixFactory(Wcols, Wrows);
  LinearOperator<double> W = Wmf->createMatrix();
  RCP<DenseSerialMatrix> Wptr = DenseSerialMatrix::getConcretePtr(W);
  vector<vector<double> > WVals = {{6, 10, 9, 7, 6, 2, 8, 9}, 
				   {6, 7, 5, 8, 2, 4, 5, 1}, 
				   {7, 4, 5, 2, 3, 2, 8, 3}, 
				   {3, 2, 2, 4, 5, 9, 1, 3}, 
				   {5, 2, 6, 3, 6, 9, 2, 2}, 
				   {6, 8, 7, 2, 1, 6, 10, 1}, 
				   {3, 1, 10, 4, 2, 4, 1, 9}, 
				   {3, 1, 9, 7, 1, 1, 9, 1}};
  Wptr->fill(WVals);

  // LinearOperator<double> Wt = Wmf->createMatrix();
  // RCP<DenseSerialMatrix> Wtptr = DenseSerialMatrix::getConcretePtr(Wt);
  // vector<vector<double> > WtVals = {{6, 6, 7, 3, 5, 6, 3, 3}, 
  // 				    {10, 7, 4, 2, 2, 8, 1, 1}, 
  // 				    {9, 5, 5, 2, 6, 7, 10, 9}, 
  // 				    {7, 8, 2, 4, 3, 2, 4, 7}, 
  // 				    {6, 2, 3, 5, 6, 1, 2, 1}, 
  // 				    {2, 4, 2, 9, 9, 6, 4, 1}, 
  // 				    {8, 5, 8, 1, 2, 10, 1, 9}, 
  // 				    {9, 1, 3, 3, 2, 1, 9, 1}};
  // Wtptr->fill(WtVals);

  // Create B = W^T S W
  //  LinearOperator<double> B = denseMatrixMatrixProduct(Wt, epetraDenseProduct(S,W));

  // Trying to check the transpose operator.
  LinearOperator<double> B = denseMatrixMatrixProduct(W.transpose(), epetraDenseProduct(S,W));
  // Start by printing out B to see it is not the 0 matrix. Compare with Mathematica

  LinearOperator<double> U2;
  LinearOperator<double> V2t;
  Vector<double> sigma2;

  // Find the SVD
  denseSVD(B, U2, sigma2, V2t);

  LinearOperator<double> negV2t = Wmf->createMatrix();
  RCP<DenseSerialMatrix> negV2tptr = DenseSerialMatrix::getConcretePtr(negV2t);
  RCP<DenseSerialMatrix> V2tptr = DenseSerialMatrix::getConcretePtr(V2t);
  RCP<DenseSerialMatrix> U2ptr = DenseSerialMatrix::getConcretePtr(U2);
  for(int i=0; i < V2t.range().dim(); i++)
    {
      for(int j=0; j < V2t.domain().dim(); j++)
	{
	  negV2tptr->setElement(i,j, -1.0*V2tptr->getElement(i,j));
	}
    }

  //cout << "Value of U - V: " << U2.addedOperator(V2t) << endl;
  cout << "Here is the matrix of IP for U*Ut: " << endl
       << denseMatrixMatrixProduct(U2,U2.transpose()) << endl;
  
  cout << "Here is the matrix of IP for U*Vt: " << endl
       << denseMatrixMatrixProduct(U2,V2t) << endl;

  cout << "Singluar values = Eigenvalues: " << endl << sigma2 << endl;

  // Go to mathematica and try to get a singular value with a few large singular matrix and a bunch of small ones, find a random orthogonal matrix U; generate a random matrix and do a QR factorization (Q will be orthognal). Form Q Sigma Q^T. I select he values for Sigma (diagonal matrix). This matrix will have the right singular values. Put that into the SVD code; check that orthogonality. the matrix I get from denseSVD may not necessarily be the same as the Q I start with. The singular values need to be the same and U and V are orthogonal for all vectors.
	
}
