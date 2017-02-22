#include "ODERHSBase.hpp"

ODERHSBase::ODERHSBase(int n)
{
    VectorType<double> vt = new SerialVectorType();
    vs_ = vt.createEvenlyPartitionedSpace(MPIComm::self(), n);
}

const VectorSpace<double>& ODERHSBase::space() const 
{return vs_;}

