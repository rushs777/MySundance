#include "KKT_Transient_Channel.hpp"

KKT_Transient_Channel::KKT_Transient_Channel(Mesh timeMesh, int Ru)
  : KKTBase(timeMesh, Ru)
{}

void KKT_Transient_Channel::stateEqn()
{

}

void KKT_Transient_Channel::adjointEqn()
{

}

void KKT_Transient_Channel::designEqn()
{

}

Expr KKT_Transient_Channel::solve()
{
  alpha_ = 10;
  return 0.0;
}
