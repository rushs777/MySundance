#ifndef SOperator_CLASS_HPP
#define SOperator_CLASS_HPP

// Sundance Includes
#include "Sundance.hpp"
//#include "PlayaDenseSerialMatrix.hpp"
//#include "VientoSnapshotIO.hpp"
//#include "denseSerialMatrixIO.hpp"

// Local Includes
#include "MathematicaConverter.hpp"


/**
 * integralOperator is a class that reflects what values the Lidar sensor returns
 * This is integration over the spatial domain only
 */
class integralOperator
{
public:
  /** Constructor */
  integralOperator(const Mesh spatialMesh, const Expr sampleDir, const QuadratureFamily quad);

  /** staticDetect returns the "true" value for the wind's velocity at node_i */
  double staticDetect(const CellFilter sensorLocation, const Expr function);

private:
  Mesh mesh_;
  Expr sampleDir_;
  QuadratureFamily quad_;
};



#endif
