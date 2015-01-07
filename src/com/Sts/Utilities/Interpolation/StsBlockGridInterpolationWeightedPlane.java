
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

/** Given dZdX and dZdY at each interpolatingPoint, compute 1/R weighted
 *  estimate of Z at point.
 */
public class StsBlockGridInterpolationWeightedPlane extends StsBlockGridInterpolation
{
    static StsBlockGridInterpolationWeightedPlane interpolator;

    static public StsBlockGridInterpolationWeightedPlane getInstance(StsModelSurface surface, StsBlockGrid blockGrid)
    {
        if(interpolator == null) interpolator = new StsBlockGridInterpolationWeightedPlane();
        interpolator.initializeGrid(surface, blockGrid);
        return interpolator;
    }

    public void initializeGrid(StsModelSurface surface, StsBlockGrid blockGrid)
    {
        if(this.surface == surface && this.blockGrid == blockGrid) return;
        if(surface.faultDistances == null && StsBlock.hasModelAuxiliarySections())
            surface.computeFaultArcLengthWeights();
        if(surface.dZdX == null) surface.computeGradients();
        super.initializeGrid(surface, blockGrid);
    }

    static public void extrapolator(StsModelSurface surface, StsBlockGrid blockGrid)
    {
        if(surface.dZdX == null) surface.computeGradients();

        if(!(interpolator instanceof StsBlockGridInterpolationWeightedPlane)) interpolator = new StsBlockGridInterpolationWeightedPlane();
        interpolator.initializeGrid(surface, blockGrid);
    }

    public float interpolate(int iCenter, int jCenter, boolean useGradient, boolean isInside, boolean debugInterpolate)
    {
        if(interpolator == null)
        {
            StsException.outputException(new StsException(StsException.FATAL, "StsInterpolateLstSqPlaneFit.interpolate() failed.",
                                            "interpolator has not been constructed."));
            return nullValue;
        }

        interpolator.initialize(iCenter, jCenter, useGradient, isInside, debugInterpolate);
        return interpolator.interpolatePoint();
    }

    float getZ(ArrayList sortedPoints)
    {
        StsBlockGridInterpolation.Point point;
        int row, col;
        float dZdX, dZdY;
        float z0, w, z, wz, wt, dx, dy;

        try
        {
            boolean debug = debugInterpolate || Main.debugPoint && iCenter == Main.debugPointRow && jCenter == Main.debugPointCol;

            int nTotalPoints = sortedPoints.size();
            if(nTotalPoints == 0) return nullValue;

            nPointsUsed = Math.min(nTotalPoints, currentMinNPoints);
            point = (StsBlockGridInterpolation.Point)sortedPoints.get(nPointsUsed-1);
            float weightMin = point.wt;
            for(int n = nPointsUsed+1; n < nTotalPoints; n++)
            {
                point = (StsBlockGridInterpolation.Point)sortedPoints.get(n);
                if(point.wt < weightMin) break;
                nPointsUsed++;
            }

            // remove the points we are not going to use
            for(int n = nTotalPoints-1; n >= nPointsUsed; n--)
                sortedPoints.remove(n);

            // useGradient is ignored

            w = 0.0f;
            wz = 0.0f;

            for(int n = 0; n < nPointsUsed; n++)
            {
                point = (StsBlockGridInterpolation.Point)sortedPoints.get(n);
                row = point.row;
                col = point.col;
                String typeName = StsParameters.getGapTypeName(blockGrid.getPointType(row, col));
                dx = -point.dx;
                dy = -point.dy;
                z0 = point.z;

                dZdX = surface.dZdX(row, col);
                dZdY = surface.dZdY(row, col);

                z = z0 + dZdX*dx + dZdY*dy;
                wt = point.wt;
                w += wt;
                wz += wt*z;

                if(debug)
                    StsException.systemDebug(this, "getZ", " point " + n + " row: " + row + " col: " + col +
                                                   " z0: " + z0 + "type: " + typeName + " wt: " + wt +
                                                   " z: " + z + " dZdX: " + dZdX + " dZdY: " + dZdY +
                                                   " dx: " + dx + " dy: " + dy);
            }
            z = wz/w;
            if(debug) StsException.systemDebug(this, "getZ", "    z: " + z);
            return z;
        }
        catch(Exception e)
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
