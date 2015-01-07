package com.Sts.Interfaces;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2002</p>
 *
 * <p>Company: </p>
 *
 * @author not attributable
 * @version 1.0
 */
public interface StsXYGridable extends StsGridable
{
    double getXOrigin();
    double getYOrigin();
    float getXSize();
    float getYSize();
    float getXMin();
    float getYMin();
    float getXInc();
    float getYInc();
    float getAngle();

    float getRowCoor(float[] xy);
    float getColCoor(float[] xy);

    float getYCoor(float rowF, float colF);
    float getXCoor(float rowF, float colF);
}
