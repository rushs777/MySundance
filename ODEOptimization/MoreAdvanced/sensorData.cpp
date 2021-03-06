#include "sensorData.hpp"

sensorData::sensorData(Array<Point> sensorLocations, Mesh spatialMesh, Mesh timeMesh, string snapshotFilenamePrefix, int nSteps, Expr sampleDir, int quadOrder)
  : allSensorLocations_(sensorLocations),
    spatialMesh_(spatialMesh),
    snapshotFilenamePrefix_(snapshotFilenamePrefix),
    nSteps_(nSteps),
    sampleDir_(sampleDir),
    quad_(new GaussianQuadrature(quadOrder))
{
  // This check ensures that the sensor locations are a subset of the spatial mesh vertices
  int numOfVertices = spatialMesh_.numCells(0);
  int Ns = allSensorLocations_.length(); // Number of sensors
  for(int i = 0; i < Ns; i++)
    {
      bool flag = false;
      for(int j = 0; j < numOfVertices; j++)
	{
	  if(allSensorLocations_[i][0] == spatialMesh_.nodePosition(j)[0] &&
	     allSensorLocations_[i][1] == spatialMesh_.nodePosition(j)[1] )
	    {
	      flag = true;
	      break;
	    }
	}
      TEUCHOS_TEST_FOR_EXCEPTION(!flag, std::runtime_error, "sensor locations are not a subset of spatial mesh vertices");	    
    }
  

}

void sensorData::create_vstar()
{
  // Create the time DiscreteSpace
  BasisFamily time_basis = new Lagrange(1);
  VectorType<double> epetraVecType = new EpetraVectorType();
  DiscreteSpace time_DS(timeMesh_, time_basis, epetraVecType);
  
  
  // This will get uForwardProblem(t_i = timeStep, x) for t_i = 0:tFinal
  Array<Expr> velocityDF;
  velocityDF.resize(nSteps_+1);
  for(int timeStep = 0; timeStep < velocityDF.length(); timeStep++)
    velocityDF[timeStep] = readSnap(snapshotFilenamePrefix_, timeStep, spatialMesh_);

  CellFilter pointFilter;
  integralOperator S(spatialMesh_,sampleDir_,quad_);

  // This takes the values from the data file and creates a time-dependent function
  // vi at each sensor location i
  // The i loop is over the sensor locations, and the j loop is over the time steps
  for(int i = 0; i < allSensorLocations_.length(); i++) // Going through each vstar[i] to make it a DiscreteFunction
    {
      Vector<double> values = time_DS.createVector();
      pointFilter = getPoint(allSensorLocations_[i]);
	  
      for(int j = 0; j < velocityDF.length(); j++) // Get the measurement values for vstar[i](t_j)
	{
	  values[j] = S.staticDetect(pointFilter, velocityDF[j]);
	  //cout << "The value of the integral eastVec*velocityDF over the point " << positionArray[i] << " at time step " << j << "  is " << values[j] << endl;	      
	}
      Expr vi = new DiscreteFunction(time_DS, values,"v["+Teuchos::toString(i)+"]");
      vstar_.append(vi);
    }
}



Expr sensorData::get_vstar()
{
  return vstar_;
}

integralOperator sensorData::get_S()
{
  integralOperator S(spatialMesh_,sampleDir_,quad_);
  return S;
}

Array<Point> sensorData::get_locations()
{
  return allSensorLocations_;
}
