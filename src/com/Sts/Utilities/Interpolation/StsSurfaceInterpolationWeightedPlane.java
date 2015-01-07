//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Utilities.Interpolation;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;

import java.util.*;

/**
 * Given dZdX and dZdY at each interpolatingPoint, compute 1/R weighted
 * estimate of Z at point.
 */
public class StsSurfaceInterpolationWeightedPlane extends StsSurfaceInterpolation
{
    static StsSurfaceInterpolationWeightedPlane interpolator;

    public StsSurfaceInterpolationWeightedPlane()
    {
    }

    static public StsSurfaceInterpolationWeightedPlane getInstance(StsSurface surface)
    {
        if (interpolator == null) interpolator = new StsSurfaceInterpolationWeightedPlane();
        interpolator.initializeGrid(surface);
        return interpolator;
    }

    protected void computeNextSpiralWeights(int i, int j, int nSpiral)
    {
    }

    public float interpolate(int iCenter, int jCenter, boolean useGradient, boolean debugInterpolate)
    {
        interpolator.initialize(iCenter, jCenter, useGradient, debugInterpolate);
        return interpolator.interpolatePoint();
    }

    protected boolean addPointOK(int i, int j)
    {
        return surface.isPointZOK(i, j);
//        return surface.isPointNotNullOrFilled(i, j);
    }

    float getZ(ArrayList sortedPoints)
    {
        Point point;
        int row, col;
        float dZdX, dZdY;
        float z0, w, z, wz, wt, dx, dy;

        try
        {
            int nTotalPoints = sortedPoints.size();
            if (nTotalPoints == 0) return nullValue;

            nPointsUsed = Math.min(nTotalPoints, currentMinNPoints);
            point = (Point) sortedPoints.get(nPointsUsed - 1);
            float weightMin = point.wt;
            for (int n = nPointsUsed + 1; n < nTotalPoints; n++)
            {
                point = (Point) sortedPoints.get(n);
                if (point.wt < weightMin) break;
                nPointsUsed++;
            }

            // remove the points we are not going to use
            for (int n = nTotalPoints - 1; n >= nPointsUsed; n--)
                sortedPoints.remove(n);

            // useGradient is ignored

            w = 0.0f;
            wz = 0.0f;

            for (int n = 0; n < nPointsUsed; n++)
            {
                point = (Point) sortedPoints.get(n);
                row = point.row;
                col = point.col;
                dx = -point.dx;
                dy = -point.dy;
                z0 = point.z;

                dZdX = surface.dZdX(row, col);
                dZdY = surface.dZdY(row, col);

                z = z0 + dZdX * dx + dZdY * dy;
                wt = point.wt;
                w += wt;
                wz += wt * z;

                if (debugInterpolate)
                {
                    String typeName = StsParameters.getGapTypeName(surface.getNullType(row, col));
                    StsMessageFiles.infoMessage(" point " + n + " row: " + row + " col: " + col + " z: " + z + "type: " + typeName + " wt: " + wt +
                        " z0: " + z0 + " dZdX: " + dZdX + " dZdY: " + dZdY + " dx: " + dx + " dy: " + dy);
                }
            }
            z = wz / w;
            if (debugInterpolate) StsMessageFiles.infoMessage("    z: " + z);
            return z;
        }
        catch (Exception e)
        {
            StsException.outputException("StsBlockGridInterpolationWeightedPlane.getZ() failed.",
                    e, StsException.WARNING);
            return nullValue;
        }
    }
/*
    protected boolean isConverged()
    {
        if(spiralHasGoodPoints) nGoodSpirals++;
        if(nGoodSpirals < minNGoodSpirals) return false;
        return weightMinimum || nPoints >= maxNPoints;
    }
*/
}