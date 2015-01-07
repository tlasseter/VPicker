package com.Sts.Rescue;
/*****************************************************************

    RescueBlockUnitPropertyGroup.java

    Rod Hanks   April 2000

******************************************************************/
public class RescueBlockUnitPropertyGroup
{
  public RescueBlockUnitPropertyGroup(String groupNameIn,
                                      RescueBlockUnit blockUnitIn)
  {
    groupNdx = CreateGroup(groupNameIn, blockUnitIn.buNdx);
  }

  private native long CreateGroup(String groupName, long buNdx);
  protected long groupNdx = 0;
}
