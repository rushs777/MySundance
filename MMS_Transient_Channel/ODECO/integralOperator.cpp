#include "integralOperator.hpp"

integralOperator::integralOperator(const Mesh mesh, const Expr sampleDir, const QuadratureFamily quad) :
  mesh_(mesh), sampleDir_(sampleDir), quad_(quad)
{}

double integralOperator::staticDetect(const CellFilter sensorLocation, const Expr function)
{
  FunctionalEvaluator pointIntegral(mesh_, Integral(sensorLocation, sampleDir_*function, quad_));
  return pointIntegral.evaluate();
}
