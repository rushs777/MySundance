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
  sensorData(Array<Point> sensorLocations, Mesh spatialMesh, Mesh timeMesh, string snapshotFilenamePrefix, int nSteps, Expr sampleDir, int quadOrder);

  void create_vstar();

  Expr get_vstar();

  integralOperator get_S();

  Array<Point> get_locations();


    
  
private:
  Array<Point> allSensorLocations_;
  Mesh spatialMesh_;
  Mesh timeMesh_;
  string snapshotFilenamePrefix_;
  int nSteps_;
  Expr sampleDir_;
  QuadratureFamily quad_;
  Expr vstar_;

};



/**
 * getPoint(Point P) returns a CellFilter that will return the cell containing
 * the P(x,y) point
 */
inline CellFilter getPoint(Point P)
{
  CellFilter vertices = new DimensionalCellFilter(0);
  CellFilter xPos = vertices.coordSubset(0, P[0]);
  CellFilter yPos = vertices.coordSubset(1, P[1]);
  CellFilter vertex = xPos.intersection(yPos);

  return vertex;
}

#endif
