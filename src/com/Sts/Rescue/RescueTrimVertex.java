package com.Sts.Rescue;
/*****************************************************************

    RescueTrimVertex.java

    Rod Hanks   April 2000

******************************************************************/

public class RescueTrimVertex extends RescuePolyLineNode
{
  public RescueTrimVertex(RescueModel model, float xIn, float yIn, float zIn)
  {
		nodeNdx = CreateTrimVertex(model.nativeNdx, xIn, yIn, zIn);
  }

	private native long CreateTrimVertex(long modelNdx, float xIn, float yIn, float zIn);
}
 
