
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

public class StsPoint2Di implements Cloneable
{
	public int x, y;

	public StsPoint2Di()
	{
	}

	public StsPoint2Di(int x, int y)
	{
		this.x = x;
        this.y = y;
	}

// CONSTRUCTOR: integer point from float values

    public StsPoint2Di(float x, float y)
	{
		this.x = (int)x;
        this.y = (int)y;
	}

// CONSTRUCTOR: integer point from float point

    public StsPoint2Di(StsPoint2D point2D)
    {
        this.x = (int)point2D.x;
        this.y = (int)point2D.y;
    }

// CONSTRUCTOR: integer point from float point + roundOff

       public StsPoint2Di(StsPoint2D point2D, float roundOff)
    {
        this.x = (int)(point2D.x + roundOff);
        this.y = (int)(point2D.y + roundOff);
    }

// METHOD: copyTo

    public void copyTo(StsPoint2Di point)
    {
        point.x = x;
        point.y = y;
    }
}
