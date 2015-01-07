
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

import javax.media.opengl.*;
import java.util.*;

public class StsGridLine implements StsEdgeLinkable
{
    int rowOrCol;
    int rowCol;
    StsSurfaceGridable grid;
    int first = -1;
    int last = -2;
    StsEdgeLoopRadialGridLink firstLink;
    StsEdgeLoopRadialGridLink lastLink;
    int lowerBound = -1; // lower bound of columns (if ROW) or rows (if COL)
    int upperBound = -2; // upper bound of columns (if ROW) or rows (if COL)

    static final int ROW = StsParameters.ROW;
    static final int COL = StsParameters.COL;

    static final int PLUS = StsParameters.PLUS;
    static final int MINUS = StsParameters.MINUS;

    public StsGridLine()
    {
    }

    public StsGridLine(StsSurfaceGridable grid, int rowOrCol, int rowCol, StsEdgeLoopRadialGridLink firstLink, StsEdgeLoopRadialGridLink lastLink)
    {
        this.grid = grid;
        this.rowOrCol = rowOrCol;
        this.rowCol = rowCol;
        this.firstLink = firstLink;
        this.lastLink = lastLink;

        this.first = StsMath.above(firstLink.getCrossingRowColF(rowOrCol));
        this.last = StsMath.below(lastLink.getCrossingRowColF(rowOrCol));

        setBounds();
    }

    private void setBounds()
    {
        float crossingRowColF;

        if(firstLink != null)
        {
            crossingRowColF = firstLink.getCrossingRowColF(rowOrCol);
            lowerBound = StsMath.floor(crossingRowColF);
        }
        else
            lowerBound = first;

        if(lastLink != null)
        {
            crossingRowColF = lastLink.getCrossingRowColF(rowOrCol);
            upperBound = StsMath.ceiling(crossingRowColF);
        }
        else
            upperBound = last;
    }

    public void setLast(int last) { this.last = last; setBounds(); }
    public int getLast() { return last; }
    public void setFirst(int first) { this.first = first; setBounds(); }
    public int getFirst() { return first; }
    public void setLastLink(StsEdgeLoopRadialGridLink lastLink) { this.lastLink = lastLink; setBounds(); }
    public StsEdgeLoopRadialGridLink getLastLink() { return lastLink; }
    public void setFirstLink(StsEdgeLoopRadialGridLink firstLink) { this.firstLink = firstLink; setBounds(); }
    public StsEdgeLoopRadialGridLink getFirstLink() { return firstLink; }
    public int getRowCol() { return rowCol; }
    public int getRowOrCol() { return rowOrCol; }
    public int getLowerBound() { return lowerBound; }
    public int getUpperBound() { return upperBound; }

    public StsGridSectionPoint getFirstGridPoint()
    {
        if(firstLink == null) return null;
        else return firstLink.getPoint();
    }

    public StsGridSectionPoint getLastGridPoint()
    {
        if(lastLink == null) return null;
        else return lastLink.getPoint();
    }

    public StsSurfaceVertex getPrevVertex()
    {
        if(firstLink == null) return null;
        return firstLink.getPoint().getVertex();
    }

    public StsSurfaceVertex getNextVertex()
    {
        if(lastLink == null) return null;
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
        int nPoints = 0;
        int row, col;
        StsPoint point;
        StsGridSectionPoint gridPoint;
        try
        {
            if(firstLink != null) nPoints++;
            if(lastLink != null) nPoints++;
            if(first <= last) nPoints += (last - first + 1);

            StsList edgePoints  = new StsList(nPoints, 10);

            int n = 0;
            if(firstLink != null) edgePoints.add(firstLink.getPoint());

            if(first <= last)
            {
                if(rowOrCol == ROW)
                {
                    row = rowCol;
                    for(col = first; col <= last; col++)
                    {
                        point = grid.getPoint(row, col);
                        gridPoint = new StsGridSectionPoint(point, row, col, null, null, false);
                        edgePoints.add(gridPoint);
                    }
                }
                else
                {
                    col = rowCol;
                    for(row = first; row <= last; row++)
                    {
                        point = grid.getPoint(row, col);
                        gridPoint = new StsGridSectionPoint(point, row, col, null, null, false);
                        edgePoints.add(gridPoint);
                    }
                }
            }

            if(lastLink != null) edgePoints.add(lastLink.getPoint());

            return edgePoints;
        }
        catch(Exception e)
        {
            StsException.outputException("StsGridLine.getEdgePoints() failed. " +
                StsParameters.rowCol(rowOrCol) + rowCol, e, StsException.WARNING);

            return null;
        }
    }

    public String getLabel()
    {
        return new String("StsGridLine-" + StsParameters.rowCol(rowOrCol) + rowCol + " on: " + grid.getLabel());
    }

    public boolean delete() { return true; }

    public float[][] getXYZPoints()
    {
        int nPoints = 0;
        int row, col;
        try
        {
            if(firstLink != null) nPoints++;
            if(lastLink != null) nPoints++;
            if(first <= last) nPoints += (last - first + 1);

            float[][] xyzPoints = new float[nPoints][];

            int n = 0;
            if(firstLink != null) xyzPoints[n++] = firstLink.getXYZ();

            if(first <= last)
            {
                if(rowOrCol == ROW)
                {
                    row = rowCol;
                    for(col = first; col <= last; col++)
                        xyzPoints[n++] = grid.getXYZorT(row, col);
                }
                else
                {
                    col = rowCol;
                    for(row = first; row <= last; row++)
                        xyzPoints[n++] = grid.getXYZorT(row, col);
                }
            }

            if(lastLink != null) xyzPoints[n++] = lastLink.getXYZ();

            return xyzPoints;
        }
        catch(Exception e)
        {
            StsException.outputException("StsGridLine.getXYZPoints() failed. " +
                StsParameters.rowCol(rowOrCol) + rowCol, e, StsException.WARNING);

            return null;
        }
    }

    public float getFirstF()
    {
        if(firstLink != null)
            return firstLink.getCrossingRowColF(rowOrCol);
        else
            return (float)first;
    }

    public float getLastF()
    {
        if(lastLink != null)
            return lastLink.getCrossingRowColF(rowOrCol);
        else
            return (float)last;
    }

    public boolean isOK()
    {
        return first >= 0 && last >= 0;
    }

    public final class GridLineComparator implements Comparator
    {
        int rowOrCol;

        GridLineComparator(int rowOrCol)
        {
            this.rowOrCol = rowOrCol;
        }

        public int compare(Object o1, Object o2)
        {
            StsGridLine p1 = (StsGridLine)o1;
            StsGridLine p2 = (StsGridLine)o2;

            int rowCol1 = p1.getRowCol();
            int rowCol2 = p2.getRowCol();
            if(rowCol1 > rowCol2) return 1;
            else if(rowCol1 < rowCol2) return -1;
            else return 0;
        }
    }

    public void draw(GL gl)
    {
        try
        {
            float[][] xyzPoints = getXYZPoints();
            int nPoints = xyzPoints.length;

            gl.glBegin(GL.GL_LINE_STRIP);

            for(int n = 0; n < nPoints; n++)
                gl.glVertex3fv(xyzPoints[n], 0);

            gl.glEnd();
        }
        catch(Exception e)
        {
            StsException.outputException("StsGridLine.draw() failed. " +
                "For: " + getLabel(), e, StsException.WARNING);
        }
    }

     public String toString()
     {
         String firstLinkString = StsToolkit.getString(firstLink);
         String lastLinkString = StsToolkit.getString(lastLink);
         return "StsGridLine from first link " + firstLinkString + " to last link " + lastLinkString + " on: " + grid.getLabel();
     }
 }








