package com.Sts.Rescue;
/*****************************************************************

    RescuePolyLineNode.java

    Rod Hanks   April 2000

******************************************************************/

public class RescuePolyLineNode
{
  public RescuePolyLineNode(float xIn, float yIn, float zIn)
  {
		nodeNdx = CreatePolyNode(xIn, yIn, zIn);
  }

	protected RescuePolyLineNode()
	{
	}

	public void SetUVValue(RescueSurface surface, float u, float v)
	{
		SetUVValuePrimitive(nodeNdx, surface.surfaceNdx, u, v);
	}

	private native long CreatePolyNode(float xIn, float yIn, float zIn);
	private native void SetUVValuePrimitive(long nodeNdx, long surfaceNdx, float u, float v);

  protected long nodeNdx = 0;
}
 
