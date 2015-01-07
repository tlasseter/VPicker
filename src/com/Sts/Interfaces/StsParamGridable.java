
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Interfaces;

import com.Sts.Types.*;

public interface StsParamGridable
{
    int getNRows();
    int getNCols();
    int getRowMin();
    int getRowMax();
    int getColMin();
    int getColMax();

    float getRowCoor(float[] xyz);
    float getColCoor(float[] xyz);

    StsPoint getPoint(int row, int col);
    public float[] getXYZ(int row, int col);

    public void checkConstructGridNormals();
    public float[] getNormal(int row, int col);
    float[] getNormal(float rowF, float colF);

    public String getLabel();
}
