
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;
import com.Sts.Types.*;

public class StsDirectedEdge
{
    protected StsEdgeLinkable edge;
    protected int direction;

    static public final int PLUS = StsParameters.PLUS;
    static public final int MINUS = StsParameters.MINUS;

    public StsDirectedEdge()
    {
    }

    public StsDirectedEdge(StsEdgeLinkable edge, int direction)
    {
        this.edge = edge;
        this.direction = direction;
    }

    public boolean initialize(StsModel model) { return true; }

    public StsEdgeLinkable getEdge() { return edge; }
    public int getDirection() { return direction; }
/*
    public void addToLinkedGrid(StsEdgeLoopRadialLinkGrid linkedGrid, boolean skipFirst)
    {
        edge.addToLinkedGrid(linkedGrid, direction, skipFirst);
    }
*/

/*
    public void addToLinkGrid(StsEdgeLoopRadialLinkGrid linkGrid, boolean skipFirst, int ID)
    {
        edge.addToLinkGrid(linkGrid, direction, skipFirst, ID);
    }
*/
    public StsList getEdgePoints()
    {
        StsList points = edge.getGridEdgePointsList();
        if( direction == PLUS ) return points;

        int len = points.getSize();
        StsObjectList reverse = new StsObjectList(len);
        for( int i=0; i<len; i++ )
            reverse.add(points.getElement(len-i-1));
        return reverse;
    }


/*
    public void setSectionRowColF(StsSection section)
    {
        StsObjectList edgePoints = edge.getEdgePoints();
        int nPnts = edgePoints.getSize();

        for(int n = 0; n < nPnts; n++)
        {
            StsGridSectionPoint edgePoint = (StsGridSectionPoint)edgePoints.getElement(n);
            edgePoint.setSectionRowOrCol(section);
        }
    }

    public void setSectionConnects(StsSection section)
    {
        section.setSectionConnects(edge, direction);
    }

    public void setGridConnects()
    {
        StsObjectList edgePoints = edge.getEdgePoints();
        int nPnts = edgePoints.getSize();

        StsGridSectionPoint point0 = null, point1 = null;

        point1 = (StsGridSectionPoint)edgePoints.getElement(0);
        for(int n = 1; n < nPnts; n++)
        {
            point0 = point1;
            point1 = (StsGridSectionPoint)edgePoints.getElement(n);
            point0.setGridConnect(point1, direction);
        }
        point1.setGridConnect(point0, -direction);
    }

    public void clearSectionConnects(StsSection section)
    {
        section.clearSectionConnects(edge);
    }
*/
    public StsGridSectionPoint getFirstPoint()
    {
        if(direction == PLUS)
            return (StsGridSectionPoint)edge.getGridEdgePointsList().getFirst();
        else if(direction == MINUS)
            return (StsGridSectionPoint)edge.getGridEdgePointsList().getLast();
        else
            return null;
    /*
        if(direction == PLUS)
            return edge.getPrevVertex().getSurfacePoint();
        else if(direction == MINUS)
            return edge.getNextVertex().getSurfacePoint();
        else
            return null;
    */
    }

    public StsGridSectionPoint getLastPoint()
    {
        if(direction == PLUS)
            return (StsGridSectionPoint)edge.getGridEdgePointsList().getLast();
        else if(direction == MINUS)
            return (StsGridSectionPoint)edge.getGridEdgePointsList().getFirst();
        else
            return null;
    /*
        if(direction == MINUS)
            return edge.getPrevVertex().getSurfacePoint();
        else if(direction == PLUS)
            return edge.getNextVertex().getSurfacePoint();
        else
            return null;
    */
    }

    public boolean delete()
    {
        edge.delete();
        return true;
    }
    public String toString()
    {
        return edge.getLabel() + "direction " + direction;
    }
}
