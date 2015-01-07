package com.Sts.Rescue;
/*****************************************************************

    RescueProperty.java

    Rod Hanks   April 2000

******************************************************************/
public class RescueProperty
{
  public RescueProperty(RescueBlockUnitPropertyGroup parentGroupIn,
                        String propertyNameIn, 
                        String propertyTypeIn, 
                        String unitOfMeasureIn,
                         float nullValue, 
                         float[][][] valueArray)
  {
/*
  Order of array indices is, quixotically, [layer][row][col], which is kij.
*/
    propNdx = CreateProperty(parentGroupIn.groupNdx,
                             propertyNameIn,
                             propertyTypeIn,
                             unitOfMeasureIn,
                             nullValue,
                             valueArray);
  }

  public RescueArray Data()
  {
    if (data == null)
    {
      long dataNdx = GetDataFor(propNdx);
      data = new RescueArray(dataNdx);
    }
    return data;
  }

  private native long CreateProperty(long groupNdx, String propertyName,
                                     String propertyType,
                                     String unitOfMeasure,
                                     float nullValue,
                                     float[][][] valueArray);
  private native long GetDataFor(long propNdx);

  private RescueArray data = null;
  protected long propNdx = 0;
}

