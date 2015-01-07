package com.Sts.Rescue;
/*****************************************************************

    RescueBlockUnit.java

    Rod Hanks   April 2000

******************************************************************/
public class RescueBlockUnit
{
  public RescueBlockUnit(int orientation,
                         RescueBlock parentBlock,
                         RescueUnit parentUnit,
                         float i_origin, float i_step,
                         int i_lowbound, int i_count,
                         float j_origin, float j_step,
                         int j_lowbound, int j_count,
                         int k_lowbound, int k_count,
                         float missingValue,
                         float rotation, RescueVertex vertex)
  {
    long vertexNdx = 0;
    if (vertex != null)
    {
      vertexNdx = vertex.vertexNdx;
    }
    buNdx = CreateBlock1(orientation, parentBlock.blockNdx, parentUnit.unitNdx,
                        i_origin, i_step, i_lowbound, i_count,
                        j_origin, j_step, j_lowbound, j_count,
                        k_lowbound, k_count, missingValue, rotation, vertexNdx);
  }

  public RescueBlockUnit(int orientation,
                  RescueBlock parentBlock,
                  RescueUnit parentUnit,
                  float i_origin, float i_step,
                  int i_lowbound, int i_count,
                  float j_origin, float j_step,
                  int j_lowbound, int j_count,
                  int k_lowbound, int k_count,
                  float missingValue,
                  RescueReferenceSurface topSurfaceIn, float topOffsetIn,
                  RescueReferenceSurface bottomSurfaceIn, float  bottomOffsetIn,
                  float rotation, RescueVertex vertex)
  {
    long vertexNdx = 0;
    if (vertex != null)
    {
      vertexNdx = vertex.vertexNdx;
    }
    buNdx = CreateBlock2(orientation, parentBlock.blockNdx, parentUnit.unitNdx,
                        i_origin, i_step, i_lowbound, i_count,
                        j_origin, j_step, j_lowbound, j_count,
                        k_lowbound, k_count, missingValue,
                        topSurfaceIn.surfaceNdx, topOffsetIn, 
                        bottomSurfaceIn.surfaceNdx, bottomOffsetIn, 
                        rotation, vertexNdx);
  }

  public RescueBlockUnit(int orientation,
                  RescueBlock parentBlock,
                  RescueUnit parentUnit,
                  float i_origin, float i_step,
                  int i_lowbound, int i_count,
                  float j_origin, float j_step,
                  int j_lowbound, int j_count,
                  int k_lowbound, int k_count,
                  float missingValue,
                  RescueReferenceSurface referenceSurfaceIn, float referenceOffsetIn, 
                  float thicknessIn, int onOffLapIn,
                  float rotation, RescueVertex vertex)
  {
    long vertexNdx = 0;
    if (vertex != null)
    {
      vertexNdx = vertex.vertexNdx;
    }
    buNdx = CreateBlock3(orientation, parentBlock.blockNdx, parentUnit.unitNdx,
                        i_origin, i_step, i_lowbound, i_count,
                        j_origin, j_step, j_lowbound, j_count,
                        k_lowbound, k_count, missingValue,
                        referenceSurfaceIn.surfaceNdx, referenceOffsetIn, 
                        thicknessIn,  onOffLapIn, 
                        rotation, vertexNdx);
  }

  public RescueMacroVolume MacroVolume()
  {
    if (macroVolume == null)
    {
      long macroNdx = MacroVolumeFor(buNdx);
      macroVolume = new RescueMacroVolume(macroNdx);
    }
    return macroVolume;
  }

  public RescueGeometry GridGeometry()
  {
    if (gridGeometry == null)
    {
      long gridGeomNdx = GridGeometryFor(buNdx);
      gridGeometry = new RescueGeometry(gridGeomNdx);
    }
    return gridGeometry;
  }

  public void SetSurfaceAboveMe(RescueBlockUnitHorizonSurface surface)
  {
    long surfNdx = 0;
    if (surface != null)
    {
      surfNdx = surface.surfaceNdx;
    }
    SetSurfaceAboveMePrimitive(buNdx, surfNdx);
  }

  public void SetSurfaceBelowMe(RescueBlockUnitHorizonSurface surface)
  {
    long surfNdx = 0;
    if (surface != null)
    {
      surfNdx = surface.surfaceNdx;
    }
    SetSurfaceBelowMePrimitive(buNdx, surfNdx);
  }


  private native long CreateBlock1(int orientation, long blockNdx, long unitNdx,
                                  float i_origin, float i_step,
                                   int i_lowbound, int i_count,
                                   float j_origin, float j_step,
                                   int j_lowbound, int j_count,
                                   int k_lowbound, int k_count,
                                   float missingValue,
                                   float rotation, long vertexNdx);

  private native long CreateBlock2(int orientation, long blockNdx, long unitNdx,
                                  float i_origin, float i_step,
                                   int i_lowbound, int i_count,
                                   float j_origin, float j_step,
                                   int j_lowbound, int j_count,
                                   int k_lowbound, int k_count,
                                   float missingValue,
                                   long topSurfaceNdx, float topOffset,
                                   long bottomSurfaceNdx, float bottomOffset,
                                   float rotation, long vertexNdx);

  private native long CreateBlock3(int orientation, long blockNdx, long unitNdx,
                                  float i_origin, float i_step,
                                   int i_lowbound, int i_count,
                                   float j_origin, float j_step,
                                   int j_lowbound, int j_count,
                                   int k_lowbound, int k_count,
                                   float missingValue,
                                   long refSurfaceNdx, float refOffset,
                                   float thickness, int onOffLap,
                                   float rotation, long vertexNdx);
  private native long MacroVolumeFor(long buNdx);
  private native long GridGeometryFor(long buNdx);
  private native void SetSurfaceAboveMePrimitive(long buNdx, long surfaceNdx);
  private native void SetSurfaceBelowMePrimitive(long buNdx, long surfaceNdx);

  private RescueMacroVolume macroVolume = null;
  private RescueGeometry gridGeometry = null;
  protected long buNdx = 0;
}
 
