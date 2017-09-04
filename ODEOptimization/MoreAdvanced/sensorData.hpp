#ifndef SENSORDATA_HPP
#define SENSORDATA_HPP

// Sundance Includes
#include "Sundance.hpp"
#include "integralOperator.hpp"
#include "VientoSnapshotIO.hpp" // for readSnap

// Local Includes
#include "MathematicaConverter.hpp"

/**
 * The class sensorData creates an expression, vstar, which is a time-dependent vector-valued 
 * function which encompasses the discrete measurement values obtained from a sensor unit(s)
 */
class sensorData
{
public:
  /** Constructor */
  sensorData(Array<Point> sensorLocations, Mesh spatialMesh, DiscreteSpace time_D, string snapshotFilenamePrefix, int nSteps, Expr sampleDir, int quadOrder);

  void create_vstar();

  Expr get_vstar();

    
  
private:
  Array<Point> allSensorLocations_;
  Mesh spatialMesh_;
  DiscreteSpace time_DS_;
  string snapshotFilenamePrefix_;
  int nSteps_;
  Expr sampleDir_;
  QuadratureFamily quad_;
  Expr vstar_;

};



#endif
