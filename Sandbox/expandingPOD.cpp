#include "Sundance.hpp"
#include "PlayaSerialVectorType.hpp"
#include "PlayaDenseSerialMatrix.hpp"
#include "PlayaDenseMatrixMatrixProduct.hpp"
#include "POD_SVD.hpp"
//#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaEpetraMatrixOps.hpp"


using std::cout;
using std::endl;
using std::vector;

int main(int argc, char *argv[])
{
  
  int problem = 1;
  Sundance::setOption("problem", problem,"Select which problem to run; 1 is a toy problem, 2 is a matrix for which the singular values > 0, 3 is a matrix for which only a few singular values are noticeably larger than 0");

  int verbosity = 1;
  Sundance::setOption("verbosity",verbosity,"Select the level of verbosity for output");
  
  Sundance::init(&argc, &argv);

  
  /*
   * Things I need: snapshot matrix W, DiscreteSpace ds for basis and mesh
   */

  // Defien the VectorTypes
  Playa::VectorType<double> epetraVecType = new Playa::EpetraVectorType();
  Playa::VectorType<double> serialVecType = new Playa::SerialVectorType();

  // Create the objects U, sigma, and Vt that denseSVD will store the results of the SVD in
  LinearOperator<double> U;
  Vector<double> sigma;
  LinearOperator<double> Vt;


  
  /*
   * Problem 1: SVD on a 3x3 matrix A = {{26,38,50},{38,56,74},{50,74,98}}
   * This matrix is formed by W^T S W where S is symmetric, positive definite
   * Result: From the SVD A = U Sigma Vt,  U == V, U U^T = I,
   * and sigma = {179.599, 0.400893, 0} >= 0
   */
  if (problem == 1)
    {
      // Create the matrix A
      VectorSpace<double> cols = serialVecType.createEvenlyPartitionedSpace(MPIComm::world(), 3);
      VectorSpace<double> rows = serialVecType.createEvenlyPartitionedSpace(MPIComm::world(), 3);
      RCP<MatrixFactory<double> > mf = serialVecType.createMatrixFactory(cols, rows);
      LinearOperator<double> A = mf->createMatrix();
      RCP<DenseSerialMatrix> Aptr = DenseSerialMatrix::getConcretePtr(A);
      vector<vector<double> > AVals = {{26,38,50},{38,56,74},{50,74,98}};
      Aptr->fill(AVals);

      // Get the SVD of A
      denseSVD(A,U,sigma,Vt);
    }
  else
    {
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

      if (problem == 2)
	{
	  /*
	   * Problem 2: Perform an SVD where the matrix S is a mass matrix generated
	   * by a LinearProblem
	   * Result: Worked correctly: U == V, U U^T = I
	   * sigma = {407.085, 21.6058, 11.0476, 4.90038, 3.45931, 1.68362, 0.490138, 0.311898}
	   */


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

	  // Find the SVD
	  denseSVD(B, U, sigma, Vt);

	  LinearOperator<double> negVt = Wmf->createMatrix();
	  RCP<DenseSerialMatrix> negVtptr = DenseSerialMatrix::getConcretePtr(negVt);
	  RCP<DenseSerialMatrix> Vtptr = DenseSerialMatrix::getConcretePtr(Vt);
	  RCP<DenseSerialMatrix> Uptr = DenseSerialMatrix::getConcretePtr(U);
	  for(int i=0; i < Vt.range().dim(); i++)
	    {
	      for(int j=0; j < Vt.domain().dim(); j++)
		{
		  negVtptr->setElement(i,j, -1.0*Vtptr->getElement(i,j));
		}
	    }
      
	}
      else if (problem == 3)
	{
	  /*
	   * Problem 3: Perform the SVD on W^T Sigma W where
	   * W is a random orthogonal matrix generated by a QR factorization
	   * and Sigma is chosen so that only the first three values are noticably
	   * larger than zero
	   * Sigma = {100, 30, 10, 1e-4, 1e-6, 1e-8, 1e-10, 0}
	   * Results: U U^T = I; however U V^T \neq I
	   * Sigma = {100, 30, 10, 1e-4, 1e-6, 1e-8, 9.99988e-11, 6.46783e-16}
	   */

	  VectorSpace<double> epetraVecSpace = epetraVecType.createEvenlyPartitionedSpace(MPIComm::self(), 8);
	  Vector<double> SigmaVec = epetraVecSpace.createMember();
	  SigmaVec[0] = 100;
	  SigmaVec[1] = 30;
	  SigmaVec[2] = 10;
	  SigmaVec[3] = 1e-4;
	  SigmaVec[4] = 1e-6;
	  SigmaVec[5] = 1e-8;
	  SigmaVec[6] = 1e-10;
	  SigmaVec[7] = 0.0;

	  LinearOperator<double> SigmaMatrix = makeEpetraDiagonalMatrix(SigmaVec);

	  cout << "Here is SigmaMatrix" << endl << SigmaMatrix << endl;


	  VectorSpace<double> cols = serialVecType.createEvenlyPartitionedSpace(MPIComm::world(), SigmaMatrix.domain().dim() );
	  VectorSpace<double> rows = serialVecType.createEvenlyPartitionedSpace(MPIComm::world(), SigmaMatrix.range().dim() );
	  RCP<MatrixFactory<double> > mf = serialVecType.createMatrixFactory(cols, rows);
	  LinearOperator<double> Q = mf->createMatrix();
	  RCP<DenseSerialMatrix> Qptr = DenseSerialMatrix::getConcretePtr(Q);
	  vector<vector<double> > QVals = {{0.5035088149780135,
		0.08391813582966891,0.08391813582966891,0.08391813582966891,0.16783627165933782, 
		0.4195906791483445,0.5874269508076824,0.4195906791483445},
	       {-0.38129822835167737,0.13113175112094488,0.47468726143781154,0.6121094655645583,
		0.19355240017851644,-0.30629667328250226,0.09338903308613418,0.31210324528785777},
	       {0.16193111726691464,0.08674645626396224,0.4775673947412303,0.3251819067893662,
		-0.3015890994653779,0.39788118504439074,-0.02203119515493632,-0.618618364320538},
	       {0.407273648138005,0.5035591198259239,-0.04862983065778002,0.06337134162698989,
		0.6267832462601162,-0.17817254321887901,-0.3061846250785729,-0.23627078409979813},
	       {-0.31077815795556535,-0.17043164962878807,0.0303488847523763,0.012254532924037805,
		0.3860533613139792,0.689613060498833,-0.45844130467730515,0.1962828574609559},
	       {-0.43369739940530094,0.778190020534339,-0.2572065925688535,-0.08222803613692412,
		-0.19298638353676478,0.24472528457435835,0.1904567249556291,-0.0014843451768843095},
	       {0.3521737920147645,0.21656395747486681,-0.025983206004941323,0.19635718782982814,
		-0.5139975061046583,-0.0017990759150029031,-0.5480390621896394,0.47465662714469337},
	       {-0.01715067325570413,0.18204085409877685,0.6851590036604096,-0.6804798163691642,
		-0.008991846178191813,-0.06790613038434357,-0.07200458435968668,0.15554608658802216}
	  };

	  Qptr->fill(QVals);


	  LinearOperator<double> C = denseMatrixMatrixProduct(Q.transpose(), epetraDenseProduct(SigmaMatrix,Q));
	  // Find the SVD
	  denseSVD(C, U, sigma, Vt);
	}
      else if (problem == 4)
	{
	  /*
	   * Problem 4: Perform an SVD using the code from POD_SVD
	   * Result: Worked correctly: U U^T = I
	   * sigma = {407.085, 21.6058, 11.0476, 4.90038, 3.45931, 1.68362, 0.490138, 0.311898}
	   */

	  VectorSpace<double> Wcols=serialVecType.createEvenlyPartitionedSpace(MPIComm::world(),8);
	  VectorSpace<double> Wrows=serialVecType.createEvenlyPartitionedSpace(MPIComm::world(),8);
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

	  POD_SVD podSVD(W, ds, verbosity);
	  podSVD.calculateSVD();
	  U = podSVD.get_singular_vectors();
	  Vt = U.transpose();
	  sigma = podSVD.get_lambda();

	  podSVD.calculateBasisFunctions();

	}
      else
	cout << "A valid problem number was not given" << endl;
    }



  cout << "U is : " << endl << U << endl;
  cout << "Vt is : " << endl << Vt << endl;
  
  cout << "U U^T: " << endl << denseMatrixMatrixProduct(U,U.transpose()) << endl;
  cout << "U V^T: " << endl << denseMatrixMatrixProduct(U,Vt) << endl;

  cout << "sigma is : " << endl << sigma << endl;

  // Go to mathematica and try to get a singular value with a few large singular matrix and a bunch of small ones, find a random orthogonal matrix U; generate a random matrix and do a QR factorization (Q will be orthognal). Form Q Sigma Q^T. I select he values for Sigma (diagonal matrix). This matrix will have the right singular values. Put that into the SVD code; check that orthogonality. the matrix I get from denseSVD may not necessarily be the same as the Q I start with. The singular values need to be the same and U and V are orthogonal for all vectors.
	
}
