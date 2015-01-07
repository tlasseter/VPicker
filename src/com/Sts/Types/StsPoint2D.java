
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

public class StsPoint2D implements Cloneable
{
	public float x, y;   // should really be private

	public StsPoint2D()
	{
	}

    // CONSTRUCTOR:
    public StsPoint2D(float x, float y)
	{
		this.x = x;
        this.y = y;
	}

    // ACCESSORS
    public void setX(float x) { this.x = x; }
    public float getX() { return x; }
    public void setY(float y) { this.y = y; }
    public float getY() { return y; }

    // METHOD: copyTo
    public void copyTo(StsPoint2D point)
    {
        point.x = x;
        point.y = y;
    }

    // METHOD: interpolateFromPoints
    public void interpolateFromPoints(StsPoint2D point0, StsPoint2D point1, float f)
    {
        this.x = point0.x + f*(point1.x - point0.x);
        this.y = point0.y + f*(point1.y - point0.y);
    }
}
