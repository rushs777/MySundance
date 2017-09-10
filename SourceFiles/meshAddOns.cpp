#include "meshAddOns.hpp"

void displayVertices(Mesh mesh_)
{

  Out::os() << "This mesh has " << mesh_.numCells(0) << " vertices " << endl;
  for(int i = 0; i < mesh_.numCells(0); i++)
    {
      Point vertex = mesh_.centroid(0,i);   
      Out::os() << "Here is vertex[" << i << "]: " << vertex << endl;
    }
}




Array<double> getBoundaries(Mesh mesh_)
{
  double xmax = -1000;
  double xmin = -1000;
  double ymin = -1000;
  double ymax = -1000;
  
  
  for(int i = 0; i < mesh_.numCells(0); i++)
    {
      Point vertex = mesh_.centroid(0,i);   

      if(i==0)
	{
	  xmax = vertex[0];
	  xmin = vertex[0];
	  ymax = vertex[1];
	  ymin = vertex[1];
	}
      else
	{
	  if(xmax < vertex[0])
	    xmax = vertex[0];

	  if(xmin > vertex[0])
	    xmin = vertex[0];

	  if(ymax < vertex[1])
	    ymax = vertex[1];

	  if(ymin > vertex[1])
	    ymin = vertex[1];
	}
    }

  // Order will be xmin, xmax, ymin,ymax
  Array<double> rtn;
  rtn.resize(mesh_.spatialDim()*2);
  rtn[0] = xmin;
  rtn[1] = xmax;
  rtn[2] = ymin;
  rtn[3] = ymax;

  return rtn;
}
