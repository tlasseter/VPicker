package com.Sts.Rescue;
/*****************************************************************

    RescueModel.java

    Rod Hanks   April 2000

******************************************************************/
public class RescueModel
{
	static
	{
		System.loadLibrary("rjni");
	}

  public RescueModel(String modelName, RescueCoordinateSystem coordinateSystemIn)
  {
		long csNdx = 0;
		if (coordinateSystemIn != null)
		{
			csNdx = coordinateSystemIn.csNdx;
		}
    nativeNdx = CreateModelFor(modelName, csNdx);
  }

	public boolean ArchiveModel(String pathName, boolean binary)
	{
		return ArchiveModelPrimitive(nativeNdx, pathName, binary);
	}

  public boolean ArchiveModel()
  {
    return ArchiveModelPrimitive(nativeNdx);
  }

  public boolean UnloadWireframe()
  {
    return UnloadWireframePrimitive(nativeNdx);
  }

  public void dispose()
  {
    DeleteNativeModel(nativeNdx);
  }

  public static native String BuildDate();
  private native long CreateModelFor(String modelName, long csNdx);
  private native boolean ArchiveModelPrimitive(long ndx, String pathName, boolean binary);
	private native boolean ArchiveModelPrimitive(long ndx);
  private native boolean UnloadWireframePrimitive(long ndx);
  private native void DeleteNativeModel(long ndx);
  protected long nativeNdx = 0;  
}
