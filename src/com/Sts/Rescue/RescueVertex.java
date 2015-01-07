package com.Sts.Rescue;
/*****************************************************************

    RescueVertex.java

    Rod Hanks   April 2000

******************************************************************/
public class RescueVertex
{
	static
	{
		System.loadLibrary("rjni");
	}

  public RescueVertex(String name, 
                      RescueCoordinateSystem existingCoordinateSystem,
                      double xIn, double yIn, double zIn)
  {
    long csNdx = 0;
    if (existingCoordinateSystem != null)
    {
      csNdx = existingCoordinateSystem.csNdx;
    }
    vertexNdx = CreateVertex(name, csNdx, xIn, yIn, zIn);
  }

  private native long CreateVertex(String name, long csNdx, 
                                   double x, double y, double z);
  protected long vertexNdx = 0;
}
