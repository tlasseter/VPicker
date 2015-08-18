package com.Sts.Interfaces;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version beta 1.0
 */

import com.Sts.Utilities.*;

public interface StsDistanceTransformFace
{
    int getNRows();
    int getNCols();
    float[][] initializeDistances();
    float[] getDistanceParameters();
    float doDistanceTransform(int row, int col, StsDistanceTransformPoint[] points);
}
