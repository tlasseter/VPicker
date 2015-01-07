
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import java.awt.*;

public class StsMousePoint
{
  	public int x, y;

	public StsMousePoint()
	{
	}

	public StsMousePoint(int x, int y)
	{
    	this.x = x;
        this.y = y;
	}

	public StsMousePoint(StsMousePoint mousePoint)
	{
    	this.x = mousePoint.x;
        this.y = mousePoint.y;
	}

    public StsMousePoint lowerLeftCoor(int[] viewport)
    {
        return new StsMousePoint(x, viewport[3] - y - 1);
    }

    public int getGLMouseY(int[] viewport)
    {
        return viewport[3] - y - 1;
    }

   public StsMousePoint lowerLeftCoor(Rectangle window)
    {
        return new StsMousePoint(x, window.height - y - 1);
    }

    public boolean near(StsMousePoint pickPoint, int pickSize)
    {
        int dx = Math.abs(x - pickPoint.x);
        int dy = Math.abs(y - pickPoint.y);

        return dx <= pickSize && dy <= pickSize;
    }

	public String toString() { return new String("x: " + x + " y: " + y); }
}
