
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.DB.*;
import com.Sts.MVC.*;

public class StsVoxelKey extends StsMainObject implements StsSerializable
{
	private float min, max;
	public StsVoxelKey()
	{
	}

    public boolean initialize(StsModel model)
    {
        return true;
    }

    // CONSTRUCTOR:
    public StsVoxelKey(float min, float max)
	{
		this.min = min;
        this.max = max;
	}

    // ACCESSORS
    public void setMin(float min) { this.min = min; }
    public float getMin() { return min; }
    public void setMax(float max) { this.max = max; }
    public float getMax() { return max; }
}
