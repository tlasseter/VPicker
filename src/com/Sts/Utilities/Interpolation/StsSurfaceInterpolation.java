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
abstract public class StsSurfaceInterpolation extends StsInterpolation
{
     /** horizon grid to be interpolated or for blockGrid interpolation, the grid to be interpolated from */
    StsSurface surface;

    static final boolean debug = false;

       /** return true if this point at ij in the source surface can be used for interpolation */
    abstract protected boolean addPointOK(int i, int j);

    public StsSurfaceInterpolation()
    {
    }

	public void initializeGrid(StsSurface surface)
	{
		if(this.surface == surface) return;
        if(surface.dZdX == null) surface.computeGradients();
        this.surface = surface;
		rowMin = 0;
		rowMax = surface.nRows-1;
		colMin = 0;
		colMax = surface.nCols-1;
        setNAvailablePoints();
    }

    boolean isPointGrid(int row, int col)
    {
        return surface.isPointNotNull(row, col);
    }

    protected void computeNextSpiralWeights(int i, int j, int nSpiral)
    {
    }

    public float interpolate(int iCenter, int jCenter, boolean useGradient, boolean debugInterpolate)
    {
        interpolator.initialize(iCenter, jCenter, useGradient, debugInterpolate);
        return interpolator.interpolatePoint();
    }

    public float interpolatePoint(int i, int j)
    {
        iCenter = i;
        jCenter = j;
        float z = nullValue;
        if(surface.isPointNull(i, j))
        {
            z = interpolatePoint();
            byte nullFillFlag = StsParameters.SURF_BOUNDARY;
            if(pointInterpolated) nullFillFlag = StsParameters.SURF_GAP_FILLED;
            if(z != nullValue) surface.setPointFilled(i, j, z, nullFillFlag);
            if(useGradient) surface.computeGradients(i, j);
        }
        return z;
    }

    void addPoint(int i, int j)
    {
        int dx, dy, distSq;
        float wt;

        boolean pointOK = addPointOK(i, j);
        if(!pointOK && !debug) return;
        float z = surface.getZorT(i, j);
        if(debug)
        {
            boolean isZNull = (z == nullValue);
            byte nullType = surface.getNullType(i, j);
            String nullTypeString = StsParameters.getSurfacePointTypeName(nullType);
            if(pointOK)
            {
                if(isZNull)
                {
                    StsException.systemDebug(this, "interpolatePoint", "point is OK: " + nullTypeString + " but z is null.");
                    return;
                }
            }
            else if(!pointOK)
            {
                if(!isZNull)
                    StsException.systemDebug(this, "interpolatePoint", "point is not OK: " + nullTypeString + " but z is not null.");
                return;
            }
        }

        dx = j - jCenter;
        dy = i - iCenter;
        distSq = dx*dx + dy*dy;
        wt = 1.0f/distSq ;
        weightSum += wt;
        weightMinimum = wt < weightMinFactor*weightSum;

        // when first point found, use it to determine minNPoints required for termination
        if(nPoints++ == 0) setMinNPoints(distSq);

        Point point = new Point(i, j, z, dx, dy, wt);
        sortedPoints.add(point);
        spiralHasGoodPoints = true;
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

    public  String getGridName() { return surface.getName(); }
/*
    protected boolean isConverged()
    {
        if(spiralHasGoodPoints) nGoodSpirals++;
        if(nGoodSpirals < minNGoodSpirals) return false;
        return weightMinimum || nPoints >= maxNPoints;
    }
*/
}