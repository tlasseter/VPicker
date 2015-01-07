
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Utilities;

public class StsPickItem
{
 	public float zMin, zMax; // scaled from 0 (nearest) to 1.0 (farthest)
    public int noNames = 0;
    public int[] names = null;

	public StsPickItem()
	{
	}
}
