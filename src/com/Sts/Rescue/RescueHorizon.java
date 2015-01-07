package com.Sts.Rescue;
/*****************************************************************

    RescueHorizon.java

    Rod Hanks   April 2000

******************************************************************/
public class RescueHorizon
{
  public RescueHorizon(String newHorizonName,
									  RescueModel parentModel)
  {
    horizonNdx = CreateHorizon(newHorizonName, parentModel.nativeNdx);
  }

  public void SetUnitAboveMe(RescueUnit existingUnit)
  {
    long unitNdx = 0;
    if (existingUnit != null)
    {
      unitNdx = existingUnit.unitNdx;
    }
    SetUnitAboveMePrimitive(horizonNdx, unitNdx);
  }

  public void SetUnitBelowMe(RescueUnit existingUnit)
  {
    long unitNdx = 0;
    if (existingUnit != null)
    {
      unitNdx = existingUnit.unitNdx;
    }
    SetUnitBelowMePrimitive(horizonNdx, unitNdx);
  }

  private native long CreateHorizon(String newHorizonName,
                                    long parentModel);
  private native void SetUnitAboveMePrimitive(long myNdx, long unitNdx);
  private native void SetUnitBelowMePrimitive(long myNdx, long unitNdx);
  protected long horizonNdx = 0;
}
 
