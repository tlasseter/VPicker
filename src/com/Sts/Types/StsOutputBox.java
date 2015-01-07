
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;

import java.io.*;

public class StsOutputBox extends StsSerialize implements Cloneable, StsSerializable, Serializable
{
    public StsPoint center = null;
    public float xDim = 0.0f;
    public float yDim = 0.0f;
    public float zDim = 0.0f;

    public StsOutputBox()  {  }
    public StsOutputBox(float x, float y, float z, float xdim, float ydim, float zdim)
    {
        center = new StsPoint(x, y, z);
        xDim = xdim;
        yDim = ydim;
        zDim = zdim;
    }

    public StsPoint getCenter()
    {
        return center;
    }

    public float getXDim()
    {
        return xDim;
    }

    public float getYDim()
    {
        return yDim;
    }

    public float getZDim()
    {
        return zDim;
    }

}




