
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.DBTypes.*;
import com.Sts.Utilities.*;

public class StsLinePoints
{
    // initialized once for all lines
    static public StsPoint2D[] linePoints;

    // input for each new line
    static int rowOrCol;
    static int rowCol;
    static boolean optimizeOK;

    // initialized for each new line
    static public int min;
    static public int max;
    static float sigma, lastSigma;

    static final int ROW = StsParameters.ROW;
    static final int COL = StsParameters.COL;

    static public void initialize(int maxPnts)
    {
        linePoints = new StsPoint2D[maxPnts];
    }

    static public void delete()
    {
        linePoints = null;
    }

    static public void initializeLine(int inRowOrCol, int inRowCol, boolean inOptimizeOK)
    {
        rowOrCol = inRowOrCol;
        rowCol = inRowCol;
        optimizeOK = inOptimizeOK;
        sigma = StsParameters.largeFloat;
        min = 0;
        max = -1;
    }

    static public boolean addPointBetterFit(int row, int col, float z)
    {
        float xy;

        if(rowOrCol == ROW)
            xy = (float)col;
        else
            xy = (float)row;

        linePoints[++max] = new StsPoint2D(xy, z);

        lastSigma = sigma;
        StsGridLineFit line = new StsGridLineFit(min, max, linePoints);
        sigma = line.lineSig;

        if(sigma >= lastSigma)
        {
             max--;
             sigma = lastSigma;
             if(optimizeOK) optimize();
             return false;
        }
        else
            return true;
    }

    /** Check if removing points at beginning of line (nearest fault) improve the
     *  line-fit; keeping removing points as long as nPnts > 3
     */
    static private void optimize()
    {
        while(max - min > 2)
        {
            min++;
            StsGridLineFit line = new StsGridLineFit(min, max, linePoints);
            sigma = line.lineSig;

            if(sigma >= lastSigma)
            {
                min--;
                sigma = lastSigma;
                break;
            }

            lastSigma = sigma;
        }
    }

    // Set points between min and max to GAP_GRID indicating they are
    // good interpolation points.
    // Set points betwen 0 and < min to GAP_NOT_FILLED, indicating they
    // have been shown in decrease lineFit sigma meaning they should be
    // interpolated.
    static public void setGoodBlockGridPoints(StsBlockGrid blockGrid)
    {
        int row, col;
        int n;

        if(rowOrCol == ROW)
        {
            row = rowCol;
            for(n = 0; n < min; n++)
            {
                col = (int)linePoints[n].getX();
                blockGrid.setPointTypeIfNot(row, col, StsParameters.GAP_NOT_FILLED, StsParameters.GAP_GRID);
            }
            for(n = min; n <= max; n++)
            {
                col = (int)linePoints[n].getX();
                // blockGrid.setPointTypeIfNot(row, col, StsParameters.GAP_GRID, StsParameters.GAP_NULL);
                blockGrid.setPointType(row, col, StsParameters.GAP_GRID);
            }
        }
        else if(rowOrCol == COL)
        {
            col = rowCol;
            for(n = 0; n < min; n++)
            {
                row = (int)linePoints[n].getX();
                blockGrid.setPointTypeIfNot(row, col, StsParameters.GAP_NOT_FILLED, StsParameters.GAP_GRID);
            }
            for(n = min; n <= max; n++)
            {
                row = (int)linePoints[n].getX();
                // blockGrid.setPointTypeIfNot(row, col, StsParameters.GAP_GRID, StsParameters.GAP_NULL);
                blockGrid.setPointType(row, col, StsParameters.GAP_GRID);
            }
        }
    }
}

