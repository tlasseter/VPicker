package com.Sts.Rescue;
/*****************************************************************

    RescueTripletArray.java

    Rod Hanks   April 2000

******************************************************************/

public class RescueTripletArray
{
  protected RescueTripletArray(long geomNdxIn)
  {
	  geomNdx = geomNdxIn;
  }

  public void AssignXValue(float[] array)
  {
	  AssignXPrimitive(geomNdx, array);
  }

  public void AssignYValue(float[] array)
  {
	  AssignYPrimitive(geomNdx, array);
  }

  public void AssignZValue(float[] array)
  {
	  AssignZPrimitive(geomNdx, array);
  }

  public void Unload()
  {
    UnloadPrimitive(geomNdx);
  }

  private native void AssignXPrimitive(long geomNdx, float[] array);
  private native void AssignYPrimitive(long geomNdx, float[] array);
  private native void AssignZPrimitive(long geomNdx, float[] array);
  private native void UnloadPrimitive(long geomNdx);

  long geomNdx = 0;
}
 
