package com.Sts.Rescue;
/*****************************************************************

    RescueMacroVolume.java

    Rod Hanks   April 2000

******************************************************************/

public class RescueMacroVolume
{
  protected RescueMacroVolume(long macroNdxIn)
  {
	  macroNdx = macroNdxIn;
  }

  public void AddBlockUnitSide(RescueBlockUnitSide toAdd)
  {
    AddBlockUnitSidePrimitive(macroNdx, toAdd.sideNdx);
  }

  public void AddInteriorSection(RescueSection toAdd)
  {
    AddInteriorSectionPrimitive(macroNdx, toAdd.surfaceNdx);
  }

  public void AddKLayerEdge(RescueEdgeSet toAdd)
  {
    AddKLayerEdgePrimitive(macroNdx, toAdd.edgeNdx);
  }

  public void SetTopEdge(RescueEdgeSet toSet)
  {
    SetTopEdgePrimitive(macroNdx, toSet.edgeNdx);
  }

  public void SetBottomEdge(RescueEdgeSet toSet)
  {
    SetBottomEdgePrimitive(macroNdx, toSet.edgeNdx);
  }

  private native void AddBlockUnitSidePrimitive(long macroNdx, long sideNdx);
  private native void AddInteriorSectionPrimitive(long macroNdx, long sectionNdx);
  private native void AddKLayerEdgePrimitive(long macroNdx, long edgeNdx);
  private native void SetTopEdgePrimitive(long macroNdx, long edgeNdx);
  private native void SetBottomEdgePrimitive(long macroNdx, long edgeNdx);
  protected long macroNdx = 0;
}
 
