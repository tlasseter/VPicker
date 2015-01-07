package com.Sts.Rescue;
/*****************************************************************

    RescueEdgeSet.java

    Rod Hanks   April 2000

******************************************************************/

public class RescueEdgeSet
{
  public RescueEdgeSet()
  {
    edgeNdx = CreateEdgeSet();
  }

  protected RescueEdgeSet(long edgeNdxIn)
  {
	  edgeNdx = edgeNdxIn;
  }

  public void AddBoundaryLoop(RescueTrimLoop toAdd)
  {
	  AddBoundary(edgeNdx, toAdd.trimNdx);
  }

  public void AddInteriorLoop(RescueTrimLoop toAdd)
  {
	  AddInterior(edgeNdx, toAdd.trimNdx);
  }

  private native void AddBoundary(long refNdx, long trimNdx);
  private native void AddInterior(long refNdx, long trimNdx);
  private native long CreateEdgeSet();

  long edgeNdx = 0;
}
 
