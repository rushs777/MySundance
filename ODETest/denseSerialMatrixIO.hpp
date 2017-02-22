#ifndef denseSerialMatrixIO_HPP
#define denseSerialMatrixIO_HPP

#include "PlayaSerialVectorType.hpp"
#include "PlayaDenseSerialMatrix.hpp"
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"


using std::vector;
using std::ofstream;
using std::ifstream;
using namespace Playa;
using namespace PlayaExprTemplates;

/** */
void writeDenseSerialMatrix(const RCP<DenseSerialMatrix>& A, string filename);

/** */
void readDenseSerialMatrix(RCP<DenseSerialMatrix>& A, string filename);


#endif
