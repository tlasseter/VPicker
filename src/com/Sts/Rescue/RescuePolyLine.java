package com.Sts.Rescue;
/*****************************************************************

    RescuePolyLine.java

    Rod Hanks   April 2000

******************************************************************/

public class RescuePolyLine
{
  public RescuePolyLine(RescueModel modelIn, RescueTrimVertex start,
                        RescueTrimVertex end)
  {
		lineNdx = CreatePolyLine(modelIn.nativeNdx, start.nodeNdx, end.nodeNdx);
  }

	public void AddPolyLineNode(RescuePolyLineNode toAdd)
	{
		AddPolyLineNodePrimitive(lineNdx, toAdd.nodeNdx);
	}

	private native long CreatePolyLine(long modelNdx, long startNdx, long endNdx);
	private native void AddPolyLineNodePrimitive(long lineNdx, long nodeNdx);

  protected long lineNdx = 0;
}
 
