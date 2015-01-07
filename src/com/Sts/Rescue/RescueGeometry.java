package com.Sts.Rescue;
/*****************************************************************

    RescueGeometry.java

    Rod Hanks   April 2000

******************************************************************/
public class RescueGeometry
{
  protected RescueGeometry(long index)
  {
    geomNdx = index;
  }

  public void AssignZValue(float[] array)
  {
    AssignZPrimitive(geomNdx, array);
  }


  public static final int R_EQUAL_AXIS = 0;						// x, y, z derived from grid origin and step or from reference surfaces.
	public static final int R_UNEQUAL_AXIS = 1;					// x and y from origin and step, z is given.
  public static final int R_COORDINATE_LINE = 2;			// straight coordinate line.
	public static final int R_COORDINATE_POLYLINE = 3;	// piece wise coordinate line.
  public static final int R_SPLIT_LINE = 4;						// split nodes (corner point or multiple coordinate line).
									   //RescueVertexType;

  public static final int R_ONLAP = 0;                // Reference surface is ABOVE the grid.
  public static final int R_OFFLAP = 1;               // Reference surface is BELOW the grid.
                      //RescueLapType;

  public static final int R_XY_ORTHOGONAL = 0;         // If VertexIs() == R_UNEQUAL_AXIS  or VertexIs() == R_EQUAL_AXIS
  public static final int R_SQUASHED_ORTHOGONAL = 1;   // If all interior vertexes are R_EQUAL_AXIS or R_UNEQUAL_AXIS but
                                                       // at least some exterior vertexes are R_COORDINATE_LINE or
                                                       // R_COORDINATE_POLYLINE.
  public static final int R_CONFORMABLE = 2;           // If neither of the above is TRUE.
                      //RescueGridType;

  public static final int R_ACTIVE = 0;                 // Cell is completely inside the BUG
  public static final int R_INACTIVE = 1;               // Cell is completely outside the BUG
  public static final int R_TRUNCATED = 2;              // Cell is partially inside and partially outside the BUG
                      //RescueCellStatus;

  protected long geomNdx = 0;
  private native void AssignZPrimitive(long geomNdx, float[] array);

}
