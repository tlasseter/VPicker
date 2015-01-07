package com.Sts.Utilities;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version beta 1.0
 */

public class StsDistanceTransformPoint
{
    public int row = -1;
    public int col = -1;
    public float distInc = 0;
    public float distance = StsParameters.largeFloat;

    static int nRows, nCols;

    static final float largeFloat = StsParameters.largeFloat;

    public StsDistanceTransformPoint(float distInc)
    {
        this.distInc = distInc;
    }

    static void setNRowsCols(int nRows_, int nCols_)
    {
        nRows = nRows_;
        nCols = nCols_;
    }

    /** returns true if not a null value or off the grid */
    float initialize(int row, int col, float[][] distances)
    {
        this.row = row;
        this.col = col;
        if(row < 0 || row >= nRows || col < 0 || col >= nCols)
            distance = largeFloat;
        else
            distance = distances[row][col];
        return distance;
    }

    public float getDistance()
    {
        if(distance == largeFloat) return distance;
        else  return distance + distInc;
    }
}