package com.Sts.Rescue;
/*****************************************************************

    RescueBlockUnitSide.java

    Rod Hanks   April 2000

******************************************************************/
public class RescueBlockUnitSide
{
  public RescueBlockUnitSide(RescueSection existingSection)
  {
    sideNdx = CreateBlockUnitSide(existingSection.surfaceNdx);
  }

  public RescueEdgeSet Edges()
  {
		if (edges == null)
		{
			long edgeNdx = GetEdgesFor(sideNdx);
			edges = new RescueEdgeSet(edgeNdx);
		}
		return edges;
  }


  private native long CreateBlockUnitSide(long sectionNdx);
  private native long GetEdgesFor(long sideNdx);
  private RescueEdgeSet edges;
  protected long sideNdx = 0;
}
 
