package com.Sts.Rescue;
/*****************************************************************

    RescueSurface.java

    Rod Hanks   April 2000

******************************************************************/
public abstract class RescueSurface
{
  public RescueSurface()
  {
  }

  public RescueEdgeSet Edges()
  {
		if (edges == null)
		{
			long edgeNdx = GetEdgesFor(surfaceNdx);
			edges = new RescueEdgeSet(edgeNdx);
		}
		return edges;
  }

  public RescueTripletArray Geometry()
  {
    if (geometry == null)
    {
      long geomNdx = GetGeometryFor(surfaceNdx);
      geometry = new RescueTripletArray(geomNdx);
    }
    return geometry;
  }

  private native long GetEdgesFor(long surfaceNdx);
  private native long GetGeometryFor(long geomNdx);

  public static int AUXILLIARY = 0;
  public static int FAULT = 1;
  public static int UNCONFORMITY = 2;
  public static int LEASE_BOUNDARY = 3;
  public static int HORIZON = 4;

  public long surfaceNdx = 0;
  private RescueEdgeSet edges = null;
  private RescueTripletArray geometry = null;
}
 
