
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.Utilities.*;

public class StsGridLinkLine implements StsEdgeLinkable
{
    private int rowOrCol;
    private int rowCol;
    private StsSurfaceGridable grid;
    StsList links;

    static final int ROW = StsParameters.ROW;
    static final int COL = StsParameters.COL;

    static final int PLUS = StsParameters.PLUS;
    static final int MINUS = StsParameters.MINUS;

    public StsGridLinkLine()
    {
    }

    public StsGridLinkLine(StsSurfaceGridable grid, int rowOrCol, int rowCol)
    {
        this.grid = grid;
        this.rowOrCol = rowOrCol;
        this.rowCol = rowCol;
        links = new StsList(10, 10);
    }

    public int getRowCol() { return rowCol; }
    public int getRowOrCol() { return rowOrCol; }

    public void addLink(StsEdgeLoopRadialGridLink link)
    {
        links.add(link);
        // link.setFlag(StsEdgeLoopRadialGridLink.USED);
    }

    public boolean isOK(int rowOrCol)
    {
        if(links == null) return false;
        int nLinks = links.getSize();
        for(int n = 0; n < nLinks; n++)
        {
            StsEdgeLoopRadialGridLink link = (StsEdgeLoopRadialGridLink)links.getElement(n);
            link.setRowColUsedFlag(rowOrCol);
        }
        return true;
    }

    public StsSurfaceVertex getPrevVertex()
    {
        StsEdgeLoopRadialGridLink firstLink = (StsEdgeLoopRadialGridLink)links.getFirst();
        return firstLink.getPoint().getVertex();
    }

    public StsSurfaceVertex getNextVertex()
    {
        StsEdgeLoopRadialGridLink lastLink = (StsEdgeLoopRadialGridLink)links.getLast();
        return lastLink.getPoint().getVertex();
    }

    public StsEdgeLinkable getNextEdge()
    {
        StsSurfaceVertex nextVertex = getNextVertex();
        if(nextVertex == null) return null;
        return nextVertex.getNextEdge();
    }

    public StsEdgeLinkable getPrevEdge()
    {
        StsSurfaceVertex prevVertex = getPrevVertex();
        if(prevVertex == null) return null;
        return prevVertex.getPrevEdge();
    }

    public StsList getEdgePointsList()
    {
        return getGridEdgePointsList();
    }

    public StsList getGridEdgePointsList()
    {
        int nPoints;

        try
        {
            nPoints = links.getSize();
            StsObjectList edgePoints  = new StsObjectList(nPoints, 10);
            for(int n = 0; n < nPoints; n++)
            {
                StsEdgeLoopRadialGridLink link = (StsEdgeLoopRadialGridLink)links.getElement(n);
                edgePoints.add(link.getPoint());
            }
            return edgePoints;
        }
        catch(Exception e)
        {
            StsException.outputException("StsGridLinkLine.getEdgePoints() failed. " +
                "for: " + getLabel(), e, StsException.WARNING);

            return null;
        }
    }

    public float[][] getXYZPoints()
    {
        int nPoints;

        try
        {
            nPoints = links.getSize();
            float[][] xyzPoints = new float[nPoints][];
            for(int n = 0; n < nPoints; n++)
            {
                StsEdgeLoopRadialGridLink link = (StsEdgeLoopRadialGridLink)links.getElement(n);
                xyzPoints[n] = link.getPoint().getXYZorT();
            }

            return xyzPoints;
        }
        catch(Exception e)
        {
            StsException.outputException("StsGridLinkLine.getXYZPoints() failed. " +
                "for: " + getLabel(), e, StsException.WARNING);

            return null;
        }
    }

    public StsGridSectionPoint getFirstGridPoint()
    {
        if(links == null) return null;
        StsEdgeLoopRadialGridLink link = (StsEdgeLoopRadialGridLink)links.getFirst();
        return link.getPoint();
    }

    public StsGridSectionPoint getLastGridPoint()
    {
        if(links == null) return null;
        StsEdgeLoopRadialGridLink link = (StsEdgeLoopRadialGridLink)links.getLast();
        return link.getPoint();
    }

    public String getLabel()
    {
        return new String("StsGridLinkLine-" + StsParameters.rowCol(rowOrCol) + rowCol +
        " on: " + grid.getLabel());
    }


    public String toString()
    {
       StsEdgeLoopRadialGridLink firstLink = (StsEdgeLoopRadialGridLink)links.getFirst();
       String firstLinkString = StsToolkit.getString(firstLink);
        StsEdgeLoopRadialGridLink lastLink = (StsEdgeLoopRadialGridLink)links.getLast();
       String lastLinkString = StsToolkit.getString(lastLink);
       return "StsGridLinkLine from first link " + firstLinkString + " to last link " + lastLinkString + " on: " + grid.getLabel();
    }
    public boolean delete() { return true; }

    public StsSurfaceGridable getGrid()
    {
        return grid;
    }
}
