
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Utilities.Interpolation;

import com.Sts.DBTypes.*;
import com.Sts.Utilities.*;

abstract public class StsBlockGridInterpolation extends StsSurfaceInterpolation
{
    StsBlockGrid blockGrid;  // block grid we are interpolating
    public boolean isInside; // point is inside of edgeLoop
    boolean hasAuxiliaryEdges = false; // edge loop includes auxiliary edges
    boolean useOutsidePoints = false; // true if point is inside and we have auxiliary edges
    boolean centerPointIsInside = true; // true if center point is inside sectionEdges
    StsModelSurface modelSurface;
    public abstract float interpolate(int iCenter, int jCenter, boolean useGradient, boolean isInside, boolean debugInterpolate);

    public void initializeGrid(StsModelSurface surface, StsBlockGrid blockGrid)
    {
        if(this.surface == surface && this.blockGrid == blockGrid) return;
        modelSurface = surface;
        if(surface.faultDistances == null && StsBlock.hasModelAuxiliarySections())
        surface.computeFaultArcLengthWeights();
        if(surface.dZdX == null) surface.computeGradients();
        this.surface = surface;
        this.blockGrid = blockGrid;
        hasAuxiliaryEdges = blockGrid.hasAuxiliaryEdges();
        rowMin = blockGrid.getRowMin();
        rowMax = blockGrid.getRowMax();
        colMin = blockGrid.getColMin();
        colMax = blockGrid.getColMax();
		hasNulls = true;
        useGradient = true;
        isInside = true;
        debugInterpolate = false;
        setNAvailablePoints();
    }

    public void initializeGrid(StsModelSurface surface, StsBlockGrid blockGrid, boolean useGradient, boolean isInside, boolean debugInterpolate)
    {
        initializeGrid(surface, blockGrid);
        this.useGradient = useGradient;
        this.isInside = isInside;
        this.debugInterpolate = debugInterpolate;
    }

    void initialize(int iCenter, int jCenter, boolean useGradient, boolean isInside, boolean debugInterpolate)
    {
        super.initialize(iCenter, jCenter, useGradient, debugInterpolate);
        this.isInside = isInside;
        useOutsidePoints = isInside && hasAuxiliaryEdges;
        try
        {
            if(useOutsidePoints) modelSurface.faultDistances[iCenter][jCenter] = 0.0f;
            if(surface != null)
                centerPointIsInside = (modelSurface.getBlockGrid(iCenter, jCenter) == blockGrid);
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "initialize", "failed at i,j: " + iCenter + " " + jCenter, e);
        }
    }

    boolean isPointGrid(int row, int col)
    {
        return blockGrid.isPointGrid(row, col);
    }

    public float interpolatePoint(int i, int j)
    {
        iCenter = i;
        jCenter = j;
        float z;
        if(debugInterpolate)
            StsException.systemDebug(this, "interpolatePoint", "row " + i + " col " + j);    
        if(surface.isPointNull(i, j))
        {
            z = interpolatePoint();
            if(z != nullValue) surface.setPointFilled(i, j, z);
        }
        else
            z = surface.getPointZ(i, j);
        return z;
    }

    protected boolean addPointOK(int i, int j)
    {
        return true;
    }

    void addPoint(int i, int j)
    {
        float faultDistance;
        int dx, dy, distSq;
        float wt;

        // Return if the surfaceGridPoint at this i,j is null.
        if(hasNulls && (blockGrid.isPointGap(i, j) || surface.isPointNull(i, j))) return;

        // useOutsidePoints is true if there are auxiliary edges and interpolated point is inside domain.
        // If point is outside domain, and there are auxiliary edges, we use only points inside the domain.

        // interpolated point is inside domain, there are auxiliary edges and this point at ij is
        // outside the domain (edgeGroup is different).

		StsBlockGrid otherBlockGrid = modelSurface.getBlockGrid(i, j);
        boolean otherPointIsInside = (otherBlockGrid == blockGrid);
        if(useOutsidePoints && centerPointIsInside && !otherPointIsInside && blockGrid.connectsBlockGrid(otherBlockGrid))
        {
            faultDistance = getFaultDistance(i, j);
            faultDistance = Math.max(0.0f, faultDistance);
            dx = j - jCenter;
            dy = i - iCenter;
            distSq = dx*dx + dy*dy;
            float faultDistSq = faultDistance*faultDistance;
            if(faultDistSq > 4*distSq) return;
            wt = 1.0f/(distSq + faultDistance);
        }
        // There are no auxiliary edges or the interpolated point is outside the domain.
        // In this case, use only points inside the block: blockGrid is the same.
        else if(otherPointIsInside && blockGrid.isPointGrid(i, j))
        {
            dx = j - jCenter;
            dy = i - iCenter;
            distSq = dx*dx + dy*dy;
            wt = 1.0f/distSq;
//            float normalWeight = normalWeights[i-rowMin][j-colMin];
//            wt = normalWeight/distSq;
        }
        else
            return;

        weightSum += wt;
        weightMinimum = wt < weightMinFactor*weightSum;

        float z = surface.getZorT(i, j);
//        float z = surface.getPointZ(i, j);

        // when first point found, use it to determine minNPoints required for termination
        if(nPoints++ == 0) setMinNPoints(distSq);

        Point point = new Point(i, j, z, dx, dy, wt);
        sortedPoints.add(point);
        spiralHasGoodPoints = true;
    }

    protected void computeNextSpiralWeights(int in, int jn, int nSpirals)
    {
        if(!useOutsidePoints || !centerPointIsInside) return;
        computeNextSpiralFaultDistances(in, jn, imin, imax, jmin, jmax, doBottom, doTop, doLeft, doRight);
    }

    float getFaultDistance(int i, int j)
    {
//        blockGrid.adjustRowColRange(i, j);
        return modelSurface.faultDistances[i][j];
    }

    void computeNextSpiralFaultDistances(int in, int jn, int iminBox, int imaxBox, int jminBox, int jmaxBox,
            boolean doBottom, boolean doTop, boolean doLeft, boolean doRight)
    {
        // set distances to largeFloat if it is outside of domain
        // and 0 if inside domain

        int imin = Math.min(iminBox, in);
        int imax = Math.max(imaxBox, in);
        int jmin = Math.min(jminBox, jn);
        int jmax = Math.max(jmaxBox, jn);

        if(doBottom)
            for(int j = jmin; j <= jmax; j++)
                initializeFaultDistance(imin, j);

        if(doTop)
            for(int j = jmin; j <= jmax; j++)
                initializeFaultDistance(imax, j);

        if(doLeft)
            for(int i = imin; i <= imax; i++)
                initializeFaultDistance(i, jmin);

        if(doRight)
            for(int i = imin; i <= imax; i++)
                initializeFaultDistance(i, jmax);

        // do bottom and top col connections
        for(int j = jmin+1; j < jmax; j++)
        {
            computeFaultDistance(imin, j, -1, 0, 1.0f);
            computeFaultDistance(imax, j, 1, 0, 1.0f);
        }

        // do left and right row connections
        for(int i = imin+1; i < imax; i++)
        {
            computeFaultDistance(i, jmin, 0, -1, 1.0f);
            computeFaultDistance(i, jmax, 0, 1, 1.0f);
        }

        // do bottom diagonal from upper left
        for(int j = jmin+1; j <= jmax; j++)
            computeFaultDistance(imin, j, -1, 1, 1.414f);

        // do bottom diagonal from upper right
        for(int j = jmin; j <= jmax-1; j++)
            computeFaultDistance(imin, j, -1, -1, 1.414f);

        // do top diagonal from lower left
        for(int j = jmin+1; j <= jmax; j++)
            computeFaultDistance(imax, j, 1, 1, 1.414f);

        // do top diagonal from lower right
        for(int j = jmin; j < jmax; j++)
            computeFaultDistance(imax, j, 1, -1, 1.414f);

        // do left diagonal from upper
        for(int i = imin; i < imax; i++)
            computeFaultDistance(i, jmin, -1, -1, 1.414f);

        // do left diagonal from lower
        for(int i = imin+1; i <= imax; i++)
            computeFaultDistance(i, jmin, 1, -1, 1.414f);

        // do right diagonal from upper
        for(int i = imin; i < imax; i++)
            computeFaultDistance(i, jmax, -1, 1, 1.414f);

        // do right diagonal from lower
        for(int i = imin+1; i <= imax; i++)
            computeFaultDistance(i, jmax, 1, 1, 1.414f);
    }

    final void initializeFaultDistance(int i, int j)
    {
        if(modelSurface.getBlockGrid(i, j) == blockGrid)
            modelSurface.faultDistances[i][j] = 0.0f;
        else
            modelSurface.faultDistances[i][j] = largeFloat;
    }

    // distance is from point at (i-di, j-dj) to (i, j)
    final void computeFaultDistance(int i, int j, int di, int dj, float dDistance)
    {
        float faultDistance, currentFaultDistance;
        float rowFaultDistance, colFaultDistance;
        float faultDistance1, faultDistance2;
        float minDistance;

        currentFaultDistance = modelSurface.faultDistances[i][j];
        if(currentFaultDistance == 0.0f) return;

        int i0 = i-di;
        int j0 = j-dj;

        float faultDistance0 = modelSurface.faultDistances[i0][j0];
        if(faultDistance0 == largeFloat) return;

        if(di == 0)
            faultDistance = getRowFaultDistance(i0, j0, dj);
        else if(dj == 0)
            faultDistance = getColFaultDistance(i0, j0, di);
        else
        {
            rowFaultDistance = getRowFaultDistance(i0, j0, dj);
            colFaultDistance = getColFaultDistance(i0, j, di);
            faultDistance1 = rowFaultDistance + colFaultDistance;

            colFaultDistance = getColFaultDistance(i0, j0, di);
            rowFaultDistance = getRowFaultDistance(i, j0, dj);
            faultDistance2 = rowFaultDistance + colFaultDistance;

            faultDistance = Math.min(faultDistance1, faultDistance2);
        }

        if(faultDistance >= largeFloat) return;

        faultDistance = faultDistance0 + faultDistance;
        modelSurface.faultDistances[i][j] = Math.min(currentFaultDistance, faultDistance);
    }

    // on row i going from j to j+dj get fault arc length distance in crossing fault
    float getRowFaultDistance(int i, int j, int dj)
    {
        float faultDistance = 0.0f;

        StsBlockGrid blockGrid0 = modelSurface.getBlockGrid(i, j);
        StsBlockGrid blockGrid1 = modelSurface.getBlockGrid(i, j+dj);
        if(blockGrid0 == null || blockGrid1 == null) return 0.0f;

        if(blockGrid0 == blockGrid1) // we didn't cross a fault
            return 0.0f;

        // we have crossed a fault: return additional distance
        if(dj > 0)
            faultDistance = getRowLinkDistance(i, j);
        else if(dj < 0)
            faultDistance = getRowLinkDistance(i, j+dj);

        if(blockGrid1 == blockGrid) // we are reentering center-point domain
            return -faultDistance;
        else
            return faultDistance;
    }

    // on col j going from i to i+di get fault arc length distance in crossing fault
    float getColFaultDistance(int i, int j, int di)
    {
        float faultDistance = 0.0f;

        StsBlockGrid blockGrid0 = modelSurface.getBlockGrid(i, j);
        StsBlockGrid blockGrid1 = modelSurface.getBlockGrid(i+di, j);
        if(blockGrid0 == null || blockGrid1 == null) return 0.0f;

        if(blockGrid0 == blockGrid1) // we didn't cross a fault
            return 0.0f;

        // we have crossed a fault: return additional distance
        if(di > 0)
            faultDistance = getColLinkDistance(i, j);
        else if(di < 0)
            faultDistance = getColLinkDistance(i+di, j);

        if(blockGrid1 == blockGrid) // we are reentering center-point domain
            return -faultDistance;
        else
            return faultDistance;
    }

    float getRowLinkDistance(int i, int j)
    {
        if(modelSurface.rowLinkDistances[i] == null) return 0.0f;
        else return modelSurface.rowLinkDistances[i][j];
    }

    float getColLinkDistance(int i, int j)
    {
        if(modelSurface.colLinkDistances[i] == null) return 0.0f;
        else return modelSurface.colLinkDistances[i][j];
    }

    public String getGridName() { return blockGrid.getLabel(); }
}
