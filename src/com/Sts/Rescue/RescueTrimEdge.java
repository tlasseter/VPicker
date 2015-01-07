package com.Sts.Rescue;
/*****************************************************************

    RescueTrimEdge.java

    Rod Hanks   April 2000

******************************************************************/

public class RescueTrimEdge
{
  public RescueTrimEdge(RescuePolyLine lineIn, int direction)
  {
		edgeNdx = CreateTrimEdge(lineIn.lineNdx, direction);
  }

  public static final int R_RIGHT_TO_LEFT = 0;
  public static final int R_LEFT_TO_RIGHT = 1;
  
	private native long CreateTrimEdge(long lineNdx, int direction);

  protected long edgeNdx = 0;
}
 
