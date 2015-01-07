package com.Sts.Rescue;
/*****************************************************************

    RescueBlockUnitHorizonSurface.java

    Rod Hanks   April 2000

******************************************************************/
public class RescueBlockUnitHorizonSurface extends RescueReferenceSurface
{
  public RescueBlockUnitHorizonSurface(int orientation,
                                        RescueHorizon parentHorizon,
                                        float i_origin, float i_step,
                                        int i_lowbound, int i_count,
                                        float j_origin, float j_step,
                                        int j_lowbound, int j_count,
                                        float missingValue,
                                        float rotation, RescueVertex vertex,
                                        int typeIn)
  {
    long vertexNdx = 0;
    if (vertex != null)
    {
      vertexNdx = vertex.vertexNdx;
    }
    surfaceNdx = CreateHorizonSurface(orientation, parentHorizon.horizonNdx, 
                                     i_origin, i_step, i_lowbound, i_count,
                                     j_origin, j_step, j_lowbound, j_count,
                                     missingValue, rotation, vertexNdx, typeIn);
  }

  public void SetBlockUnitAboveMe(RescueBlockUnit toSet)
  {
    long buNdx = 0;
    if (toSet != null)
    {
      buNdx = toSet.buNdx;
    }
    SetBlockUnitAboveMePrimitive(surfaceNdx, buNdx);
  }

  public void SetBlockUnitBelowMe(RescueBlockUnit toSet)
  {
    long buNdx = 0;
    if (toSet != null)
    {
      buNdx = toSet.buNdx;
    }
    SetBlockUnitBelowMePrimitive(surfaceNdx, buNdx);
  }

  private native long CreateHorizonSurface(int orientation,
                                            long parentHorizonNdx,
                                            float i_origin, float i_step,
                                            int i_lowbound, int i_count,
                                            float j_origin, float j_step,
                                            int j_lowbound, int j_count,
                                            float missingValue,
                                            float rotation, long vertexNdx,
                                            int typeIn);
  private native void SetBlockUnitAboveMePrimitive(long myNdx, long buNdx);
  private native void SetBlockUnitBelowMePrimitive(long myNdx, long buNdx);
}
 
