package com.Sts.Rescue;
/*****************************************************************

    RescueSection.java

    Rod Hanks   April 2000

******************************************************************/
public class RescueSection extends RescueSurface
{
  public RescueSection(int orientation,
											 String newSectionName,
											 RescueModel parentModel,
											 int typeIn,
											 int i_lowbound, int i_count,
											 int j_lowbound, int j_count,
											 float missingValue)
  {
    surfaceNdx = CreateSection(orientation, newSectionName, parentModel.nativeNdx,
                               typeIn, i_lowbound, i_count, j_lowbound, j_count,
                               missingValue);
  }
  private native long CreateSection(int orientation,
                                    String newSectionName,
                                    long parentModel,
                                    int typeIn,
                                    int i_lowbound, int i_count,
                                    int j_lowbound, int j_count,
                                    float missingValue);
}
 
