package com.Sts.Rescue;
/*****************************************************************

    RescueArray.java

    Rod Hanks   April 2000

******************************************************************/

public class RescueArray
{
  protected RescueArray(long dataNdxIn)
  {
	  dataNdx = dataNdxIn;
  }

  public void Unload()
  {
    UnloadPrimitive(dataNdx);
  }

  private native void UnloadPrimitive(long dataNdx);

  long dataNdx = 0;
}
 
