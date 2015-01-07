
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Interfaces;

public interface StsGridLinkable
{
    String getLabel();
    int getRowOrCol();
    int getRowCol();
    float getRowF();
    float getColF();
    float getRowF(StsSurfaceGridable grid);
    float getColF(StsSurfaceGridable grid);
    boolean delete();
    float[] getXYZ();
    int getID();
}