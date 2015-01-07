package com.Sts.Rescue;
/*****************************************************************

    RescueUnit.java

    Rod Hanks   April 2000

******************************************************************/
public class RescueUnit
{
  public RescueUnit(String newUnitName,
									  RescueModel parentModel)
  {
    unitNdx = CreateUnit(newUnitName, parentModel.nativeNdx);
  }

  public void SetHorizonAboveMe(RescueHorizon existingHorizon)
  {
    long horizonNdx = 0;
    if (existingHorizon != null)
    {
      horizonNdx = existingHorizon.horizonNdx;
    }
    SetHorizonAboveMePrimitive(unitNdx, horizonNdx);
  }

  public void SetHorizonBelowMe(RescueHorizon existingHorizon)
  {
    long horizonNdx = 0;
    if (existingHorizon != null)
    {
      horizonNdx = existingHorizon.horizonNdx;
    }
    SetHorizonBelowMePrimitive(unitNdx, horizonNdx);
  }

  private native long CreateUnit(String newUnitName,
                                 long parentModel);
  private native void SetHorizonAboveMePrimitive(long myNdx, long horizonNdx);
  private native void SetHorizonBelowMePrimitive(long myNdx, long horizonNdx);
  protected long unitNdx = 0;
}
 
