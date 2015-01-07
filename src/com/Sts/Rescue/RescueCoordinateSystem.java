package com.Sts.Rescue;
/*****************************************************************

    RescueCoordinateSystem.java

    Rod Hanks   April 2000

******************************************************************/
public class RescueCoordinateSystem
{
	static
	{
		System.loadLibrary("rjni");
	}

  public static final int LUF = 0;
  public static final int LUB = 1;
  public static final int LDF = 2;
  public static final int LDB = 3;
  public static final int RUF = 4;
  public static final int RUB = 5;
  public static final int RDF = 6;
  public static final int RDB = 7;
/*
  The displayOrientation shows which corner the origin is located on:

    Key:	Left/Right
	        Up/Down
	        Front/Back
*/
  public RescueCoordinateSystem(String name,
                                int displayOrientationIn,
                                RescueVertex vertexIn,
                                String Xproperty,
                                String Xuom,
                                String Yproperty,
                                String Yuom,
                                String Zproperty,
                                String Zuom)
  {

    System.loadLibrary("rjni");

    String libraryName = System.mapLibraryName("rjni");
    long vertexNdx = 0;
    if (vertexIn != null)
    {
      vertexNdx = vertexIn.vertexNdx;
    }
    csNdx = CreateCSFor(name, displayOrientationIn, vertexNdx,
                        Xproperty, Xuom, Yproperty, Yuom, Zproperty, Zuom);
  }

  private native long CreateCSFor(String name,
                                  int displayOrientationIn,
                                  long vertexNdx,
                                  String Xproperty,
                                  String Xuom,
                                  String Yproperty,
                                  String Yuom,
                                  String Zproperty,
                                  String Zuom);
  public long csNdx = 0;
}
