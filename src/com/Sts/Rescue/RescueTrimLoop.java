package com.Sts.Rescue;
/*****************************************************************

    RescueTrimLoop.java

    Rod Hanks   April 2000

******************************************************************/

public class RescueTrimLoop
{
  public RescueTrimLoop()
  {
		trimNdx = CreateTrimLoop();
  }

  public void AddLoopEdge(RescueTrimEdge toAdd)
  {
	  AddLoopEdgePrimitive(trimNdx, toAdd.edgeNdx);
  }

  private native void AddLoopEdgePrimitive(long trimNdx, long edgeNdx);
	private native long CreateTrimLoop();

  protected long trimNdx = 0;
}
 
