
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Interfaces;

import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;

public interface StsXYSurfaceGridable extends StsSurfaceGridable, StsXYGridable
{
    float getZMin();
    float getZMax();
    float getZInc();
    float[][] getPointsZ();
//    float[][] getAdjPointsZ();

	float interpolateBilinearZ(StsPoint point, boolean computeIfNull, boolean setPoint);
    float interpolateBilinearZ(StsGridPoint gridPoint, boolean computeIfNull, boolean setPoint);

    public float getComputePointZ(int row, int col);

    public boolean toggleSurfacePickingOn();
    public void toggleSurfacePickingOff();
    public String getName();
    public StsGridPoint getSurfacePosition(StsMouse mouse, boolean display, StsGLPanel3d glPanel3d);
    public void setIsVisible(boolean isVisible);
    public boolean getIsVisible();
}
