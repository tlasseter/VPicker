
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;

public class StsEdge extends StsMainObject implements StsEdgeLinkable, StsInstance3dDisplayable
{
    protected StsSection section = null;
    protected StsXYSurfaceGridable surface = null;
	/** persistent StsGridSectionPoint points along section and surface */
    protected StsObjectRefList edgePoints;
	/** First vertex on edge */
    protected StsSurfaceVertex prevVertex = null;
	/** Last vertex on edge */
    protected StsSurfaceVertex nextVertex = null;

    /** transient edgePoints used by StsSurfaceEdge only instead of persistent edgePoints. */
    transient StsObjectList tempEdgePoints = null;
    /** temporary StsGridSectionPoints along section and surface */
	transient protected StsList gridEdgePoints = null;
    //TODO remove this and directly use tempEdgePoints
    /** points extracted from StsGridSectionPoints: for drawing only */
	transient protected StsPoint[] points = null;
	/** if true, draw each points as colored square */
	transient protected boolean drawPoints = false;

    transient static StsEdge currentEdge = null;

    // Convenience copies of flags

    static final int MINUS = StsParameters.MINUS;
    static final int PLUS = StsParameters.PLUS;
    static final int NONE = StsParameters.NONE;
    static final int PLUS_AND_MINUS = StsParameters.PLUS_AND_MINUS;

    static final int X = StsPoint.X;
    static final int Y = StsPoint.Y;

    static final float LARGE_FLOAT = StsParameters.largeFloat;
    static final float nullZValue = StsParameters.nullValue;

    public static final int ROW = StsParameters.ROW;
    public static final int COL = StsParameters.COL;
    public static final int ROWCOL = StsParameters.ROWCOL;

    public static final int FIRST = StsParameters.FIRST;
    public static final int LAST = StsParameters.LAST;

    public static final int OK = StsParameters.OK;
    public static final int NOT_OK = StsParameters.NOT_OK;

// Constructors

	public StsEdge()
	{
    }

    public StsEdge(boolean persistent)
    {
        super(persistent);
    }


    public StsEdge(StsSurfaceVertex v0, StsSurfaceVertex v1, int direction)
	{
    	this();
   		v0.addEdge(this);
        v1.addEdge(this);

        if(direction > 0)
        {
		    prevVertex = v0;
		    nextVertex = v1;
        }
        else
        {
		    prevVertex = v1;
		    nextVertex = v0;
        }
	}

// Accessors

    public void setSection(StsSection section) { this.section = section; }
    public StsSection getSection() { return section; }
	public void setSurface(StsModelSurface surface) { this.surface = surface; }
    public StsXYSurfaceGridable getSurface() { return surface; }
    public StsSurfaceVertex getPrevVertex() { return prevVertex; }
    public StsSurfaceVertex getNextVertex() { return nextVertex; }
    // public void setEdgePoints(StsObjectRefList edgePoints) { this.edgePoints = edgePoints; }
    public void setDrawPoints(boolean draw) { drawPoints = draw; }
    public boolean getDrawPoints() { return drawPoints; }
    public void setPoints(StsPoint[] points) {this.points = points; }

//    public StsGrid getSurfaceGrid() { return surface.getGrid(); }
//	public StsXYSurfaceGridable getEdgeGrid() { return surface.getGrid(); }

    public StsList getGridEdgePointsList()
     {
         if(gridEdgePoints != null) return gridEdgePoints;
         return getEdgePointsList();
     }

    public StsObjectList getPersistentEdgePoints()
    {
		if(edgePoints == null) return null;
		return edgePoints.getList();
	}

    public StsPoint[] getPoints()
	{
		return points;
	}

    // included for interface compatibility
    public int getRowOrCol() { return -1; }
    public int getRowCol() { return -1; }

	public boolean setPrevVertex(StsSurfaceVertex prevVertex)
	{
		if(prevVertex == null)
		{
			StsException.systemError("StsEdge.addPrevVertex() failed: vertex is null.");
			return false;
		}
		if(this.prevVertex != null)
		{
			StsException.systemError("StsEdge.addPrevVertex() failed." +
				" prevVertex already exists: " + prevVertex.getLabel());
			return false;
		}

		this.prevVertex = prevVertex;
        prevVertex.addEdge(this);
		return true;
	}

	// adds prevVertex and edgePoint to list
	public boolean addPrevVertex(StsSurfaceVertex prevVertex)
	{
	    if(!setPrevVertex(prevVertex)) return false;
        StsObjectList pointsList = getEdgePointsList();
        if(pointsList != null && pointsList.getSize() > 0)
		{
			StsException.systemError("StsEdge.addPrevVertex() failed." +
				" edgePoints already exist: cannot add first one.");
			return false;
		}
		addEdgePoint(prevVertex.getSurfacePoint());
		return true;
	}

	public boolean setNextVertex(StsSurfaceVertex nextVertex)
	{
		if(nextVertex == null)
		{
			StsException.systemError("StsEdge.addNextVertex() failed: vertex is null");
			return false;
		}
		if(this.nextVertex != null)
		{
			StsException.systemError("StsEdge.addNextVertex() failed." +
				" nextVertex already exists: " + nextVertex.getLabel());
			return false;
		}

		this.nextVertex = nextVertex;
        nextVertex.addEdge(this);
		return true;
	}

	// adds nextVertex and edgePoint to list
	public boolean addNextVertex(StsSurfaceVertex nextVertex)
	{
	    if(!setNextVertex(nextVertex)) return false;
		addEdgePoint(nextVertex.getSurfacePoint());
		return true;
	}
/*
    public void setNextVertex(StsSurfaceVertex vertex)
    {
		if(vertex == null) return;
        nextVertex = vertex;
        vertex.addEdge(this);
		if(edgePoints == null) edgePoints = StsObjectRefList.constructor(2, 10);
		addEdgePoint(vertex.getSurfacePoint());
    }

    public void setPrevVertex(StsSurfaceVertex vertex)
    {
		if(vertex == null) return;
        prevVertex = vertex;
        vertex.addEdge(this);
		StsGridSectionPoint edgePoint = vertex.getSurfacePoint();
		if(edgePoints == null)
		    addEdgePoint(edgePoint);
		else
		    edgePoints.insertBefore(0, edgePoint);
    }
*/
    public StsEdgeLinkable getPrevEdge()
    {
        if(prevVertex == null) return null;
        return prevVertex.getOtherEdge(this);
    }
/*
    public StsEdgeLinkable getPrevEdge()
    {
        if(prevVertex == null) return null;
        return prevVertex.getPrevEdge();
    }
*/
    public StsEdgeLinkable getNextEdge()
    {
        if(nextVertex == null) return null;
        return nextVertex.getOtherEdge(this);
    }
/*
    public StsEdgeLinkable getNextEdge()
    {
        if(nextVertex == null) return null;
        return nextVertex.getNextEdge();
    }
*/
    public boolean isCurved()
	{
		StsList edgePoints = getGridEdgePointsList();
		return edgePoints != null && edgePoints.getSize() > 2;
	}


    public boolean isBoundary()
    {
        return section.isBoundary();
    }

    public boolean isFaulted()
    {
        return section.isFault();
    }

    public boolean isAuxiliary()
    {
        return section.isAuxiliary();
    }

    protected float[] getColRange()
    {
        float c0 = getPrevVertex().getSectionColF(section);
        float c1 = getNextVertex().getSectionColF(section);

        if(c0 < c1) return new float[] { c0, c1 };
        else return new float[] { c1, c0 };
    }

    public boolean initialize(StsModel model)
    {
        prevVertex.addEdge(this);
        nextVertex.addEdge(this);
        return true;
    }

    public boolean delete()
    {
        if(prevVertex != null) prevVertex.deleteEdgeFromVertex(this);
        if(nextVertex != null) nextVertex.deleteEdgeFromVertex(this);
        if(edgePoints != null) edgePoints.deleteAll();
        super.delete();
        return true;
    }

    public void deleteVertices()
    {
        if(prevVertex != null)
        {
            prevVertex.deleteEdgeFromVertex(this);
            prevVertex = null;
        }
        if(nextVertex != null)
        {
            nextVertex.deleteEdgeFromVertex(this);
            nextVertex = null;
        }
    }
/*
    public void convertEdgePointsToPersistent()
    {
        tempEdgePoints.trimToSize();
        edgePoints = tempEdgePoints.convertListToRefList(currentModel, "edgePoints", this);
		tempEdgePoints = null;
    }
*/
    /** draw a segmented dotted-line between points */

    public void display(StsGLPanel3d glPanel3d)
    {
    	StsColor stsColor;

        StsPoint[] points = getPoints();
        if(points == null) return;

        GL gl = glPanel3d.getGL();
        if(gl == null) return;

        boolean isDepth = (currentModel.getProject().getZDomain() == StsProject.TD_DEPTH);
        int zIndex;
        if(isDepth)
            zIndex = 2;
        else
            zIndex = 4;

        stsColor = StsColor.GREEN;

		glPanel3d.setViewShift(gl, StsGraphicParameters.edgeShift);
        StsGLDraw.drawDottedLine(gl, stsColor, StsGraphicParameters.edgeLineWidthHighlighted, points, zIndex);
		glPanel3d.resetViewShift(gl);

		if(drawPoints)
			StsGLDraw.drawPoints(points, StsColor.WHITE, glPanel3d,
                                 4, StsGraphicParameters.vertexShift, zIndex);
	}

    public void pick(GL gl, StsGLPanel glPanel)
    {
        float edgeWidth;

        StsPoint[] points = getPoints();
        if(points == null) return;

        StsColor color = StsColor.BLACK;

        if(section != null && StsSection.getCurrentSection() == section)
            edgeWidth = StsGraphicParameters.edgeLineWidthHighlighted;
        else
            edgeWidth = StsGraphicParameters.edgeLineWidth;

        StsGLDraw.pickEdge(gl, color, edgeWidth, points);
    }

	public StsPoint[] getPointsFromEdgePoints()
    {
        return getPointsFromEdgePoints(isDepth);
    }

	public StsPoint[] getPointsFromEdgePoints(boolean isDepth)
    {
        StsList edgePoints = getEdgePointsList();
        if(edgePoints == null) return null;
        int nEdgePoints = edgePoints.getSize();
        StsPoint[] points = new StsPoint[nEdgePoints];

        for(int i = 0; i < nEdgePoints; i++)
        {
            StsGridSectionPoint edgePoint = (StsGridSectionPoint)edgePoints.getElement(i);
            points[i] = edgePoint.getPointXYZorT(isDepth);
        }
        return points;
    }

    public void setPointsFromEdgePoints()
    {
        points = getPointsFromEdgePoints(isDepth);
    }

    public StsPoint[] getPointsFromGridEdgePoints()
    {
        if(gridEdgePoints == null) return null;
        int nEdgePoints = gridEdgePoints.getSize();
        StsPoint[] points = new StsPoint[nEdgePoints];

        for(int i = 0; i < nEdgePoints; i++)
        {
            StsGridSectionPoint gridEdgePoint = (StsGridSectionPoint)gridEdgePoints.getElement(i);
            points[i] = gridEdgePoint.getPointXYZorT();
        }
        return points;
    }

    public void setPointsFromGridEdgePoints()
    {
        points = getPointsFromGridEdgePoints();
    }
    
    public float[][] getXYZPoints()
    {
        StsObjectList edgePointsList = getEdgePointsList();
        int nPoints = edgePointsList.getSize();
        float[][] xyzPoints = new float[nPoints][];
        for(int n = 0; n < nPoints; n++)
        {
            StsGridSectionPoint edgePoint = (StsGridSectionPoint)edgePointsList.getElement(n);
            xyzPoints[n] = edgePoint.getXYZorT(isDepth);
        }
        return xyzPoints;
    }
/*
	public void addGridCrossingEdgePoints()
	{
		gridEdgePoints = section.addGridCrossingEdgePoints(this, edgePoints.getList(), null);
        points = getPointsFromGridEdgePoints();
    }
*/
    /** Using the picked points, construct interpolated points at each grid crossing.  These points are used
     *  for display (saved in edge.points) and in constructing edgeLoops on trimmed surfaces.
     * @param adjustToGrid indicates we want to adjust the interpolated points vertically to grid
     */
    public void addGridCrossingEdgePoints(boolean adjustToGrid)
	{
        if(surface != null && surface.getNRows() == 0) return;
        if(!adjustToGrid)
			gridEdgePoints = section.addGridCrossingEdgePoints(this, getEdgePointsList(), null);
        else
            gridEdgePoints = section.addGridCrossingEdgePoints(this, getEdgePointsList(), surface);

        points = getPointsFromGridEdgePoints();
    }

    public void constructSurfaceEdgePoints(StsXYSurfaceGridable grid, boolean persistent)
	{
        gridEdgePoints = section.constructSurfaceEdgePoints(this, grid, persistent);
        points = getPointsFromGridEdgePoints();
//		if(persistent) convertEdgePointsToPersistent();
	}
/*
	public void copyEdgePoints(StsEdge edge)
	{
		StsObjectList edgePoints;
        StsGridSectionPoint edgePoint0, edgePoint1, edgePoint;
        float colf0, colf1, colF;
		int nFirst, nLast, n0, n1, inc, nPoints;

		try
		{
			edgePoints = edge.getEdgePoints();
			if(edgePoints == null) return;
			int nEdgePoints = edgePoints.getSize();

			edgePoint0 = prevVertex.getSurfacePoint();
			edgePoint1 = nextVertex.getSurfacePoint();
//            edgePoint0 = (StsGridSectionPoint)edgePoints.getFirst();
//            edgePoint1 = (StsGridSectionPoint)edgePoints.getLast();

			nFirst = getEdgePointsInterval(edgePoints, edgePoint0);
			if(nFirst == -1) return;
			nLast =  getEdgePointsInterval(edgePoints, edgePoint1);
			if(nLast == -1) return;

            if(nLast >= nFirst)
            {
				n0 = nFirst+1;
				n1 = nLast;
				inc = 1;
				nPoints = n1-n0+1;
            }
			else
			{
				n0 = nFirst;
				n1 = nLast+1;
				inc = -1;
				nPoints = n0-n1+1;
			}

			tempEdgePoints = new StsObjectList(nPoints+2);
			tempEdgePoints.add(edgePoint0);
			for(int n = n0, i = 0; i < nPoints; i++, n += inc)
			{
					edgePoint = (StsGridSectionPoint)edgePoints.getElement(n);
				    tempEdgePoints.add(edgePoint);
            }
			tempEdgePoints.add(edgePoint1);
		}
		catch(Exception e)
		{
			StsException.outputException("StsEdge.copyVerticalEdgePoints() failed.",
				e, StsException.WARNING);
		}
	}
*/
	private int getEdgePointsInterval(StsList edgePoints, StsGridSectionPoint edgePoint)
	{
	    int position = edgePoint.getSectionPosition(section);
		if(position == MINUS)
			return 0;
		else if(position == PLUS)
			return edgePoints.getSize()-2;
		else
		{
		    if(section.isVertical())
				return getVerticalEdgePointsInterval(edgePoints, edgePoint);
		    else
				return getCurvedEdgePointsInterval(edgePoints, edgePoint);
		}
	}

	public int getVerticalEdgePointsInterval(StsList edgePoints, StsGridSectionPoint endPoint)
	{
		StsGridSectionPoint edgePoint;
        float colF0, colF1, colF;

		try
		{
			colF = endPoint.getColF(section);

		    int nEdgePoints = edgePoints.getSize();
			edgePoint = (StsGridSectionPoint)edgePoints.getElement(0);
			colF1 = edgePoint.getColF(section);
			for(int n = 1; n < nEdgePoints; n++)
			{
				colF0 = colF1;
				edgePoint = (StsGridSectionPoint)edgePoints.getElement(n);
				colF1 = edgePoint.getColF(section);
				if(StsMath.betweenInclusive(colF, colF0, colF1)) return n-1;
//				if(colF > colF0 && colF < colF1) return n-1;
			}
			StsException.systemError("StsEdge.getVerticalEdgePointsInterval() failed." +
				" Didn't find edgePoint on edge: " + getLabel() +
				" at sectionCol:" + colF);
			return -1;
		}
		catch(Exception e)
		{
			StsException.outputException("StsEdge.getVerticalEdgePointsInterval() failed.",
				e, StsException.WARNING);
			return -1;
		}
	}

	private int getCurvedEdgePointsInterval(StsList edgePoints, StsGridSectionPoint endPoint)
	{
		float endSectionColF, sectionColF0, sectionColF1;
		StsGridRowCol endGridRowCol, gridRowCol0, gridRowCol1;
		int nNext, inc;
		StsGridSectionPoint edgePoint;

		try
		{
			if(edgePoints == null) return -1;
			int nEdgePoints = edgePoints.getSize();
			if(nEdgePoints < 2) return -1;

			endSectionColF = endPoint.getColF(section);
			endGridRowCol = endPoint.getGridRowCol();

			edgePoint = (StsGridSectionPoint)edgePoints.getElement(0);
			sectionColF1 = edgePoint.getColF(section);
			gridRowCol1 = edgePoint.getGridRowCol();
			for(int n = 1; n < nEdgePoints; n++)
			{
				sectionColF0 = sectionColF1;
				gridRowCol0 = gridRowCol1;
				edgePoint = (StsGridSectionPoint)edgePoints.getElement(n);
				sectionColF1 = edgePoint.getColF(section);
				gridRowCol1 = edgePoint.getGridRowCol();
				if(StsMath.betweenInclusive(endSectionColF, sectionColF0, sectionColF1)) return n-1;
//				if(!StsMath.between(endSectionColF, sectionColF0, sectionColF1)) continue;
//				if(endGridRowCol.between(gridRowCol0, gridRowCol1, 0.25f)) return n-1;
			}

			StsException.systemError("StsEdge.getCurvedEdgePointsInterval() failed." +
				" Didn't find edgePoint on edge: " + getLabel() +
				" at gridRowCol:" + endGridRowCol);
			return -1;
		}
		catch(Exception e)
		{
			StsException.outputException("StsEdge.getCurvedEdgePointsInterval() failed.",
				e, StsException.WARNING);
			return -1;
		}
	}

    public StsGridSectionPoint getFirstEdgePoint()
    {
        StsObjectList edgePointsList = getEdgePointsList();
        if(edgePointsList == null) return null;
        return (StsGridSectionPoint)edgePointsList.getFirst();
    }

    public StsGridSectionPoint getLastEdgePoint()
    {
        StsObjectList edgePointsList = getEdgePointsList();
        if(edgePointsList == null) return null;
        return (StsGridSectionPoint)edgePointsList.getLast();
    }

    public StsGridSectionPoint getFirstGridPoint()
    {
        return getFirstEdgePoint();
    }

    public StsGridSectionPoint getLastGridPoint()
    {
        return getLastEdgePoint();
    }

	public StsPoint getFirstPoint()
	{
		StsObjectList edgePointsList = getEdgePointsList();
        if(edgePointsList == null) return null;
        return ((StsGridSectionPoint)edgePointsList.getFirst()).getPoint();
	}

	public StsPoint getLastPoint()
	{
		StsObjectList edgePointsList = getEdgePointsList();
        if(edgePointsList == null) return null;
        return ((StsGridSectionPoint)edgePointsList.getLast()).getPoint();
	}

	public StsPoint getSecondPoint()
	{
		StsObjectList edgePointsList = getEdgePointsList();
        if(edgePointsList == null || edgePointsList.getSize() < 2) return null;
		return ((StsGridSectionPoint)edgePointsList.getElement(1)).getPoint();
	}

	public StsPoint getSecondToLastPoint()
	{
        StsObjectList edgePointsList = getEdgePointsList();
        if(edgePointsList == null) return null;
		int nPoints = edgePointsList.getSize();
		if(nPoints < 2) return null;
		return ((StsGridSectionPoint)edgePointsList.getElement(nPoints-2)).getPoint();
	}

    public StsEdgeLoopRadialGridLink getFirstRowOrColCrossingLink(StsEdgeLoop edgeLoop)
    {
        Iterator<StsEdgeLoopRadialGridLink> edgeLinksIterator = edgeLoop.getEdgeLinksIterator(this, PLUS);
        while(edgeLinksIterator.hasNext())
        {
            StsEdgeLoopRadialGridLink link = edgeLinksIterator.next();
            if(link == null) return null;
            if(link.isRowOrCol()) return link;
        }
        return null;
    }

    public int getRowColDirection()
    {
        if(!section.isRowColAligned()) return NONE;

        StsGridSectionPoint prevEdgePoint = prevVertex.getSurfacePoint();
        StsGridSectionPoint nextEdgePoint = nextVertex.getSurfacePoint();

        if((int)prevEdgePoint.getRowF(null) == (int)nextEdgePoint.getRowF(null))
        {
            if(prevEdgePoint.getColF(null) < nextEdgePoint.getColF(null))
                return PLUS;
            else
                return MINUS;
        }
        else if((int)prevEdgePoint.getColF(null) == (int)nextEdgePoint.getColF(null))
        {
            if(prevEdgePoint.getRowF(null) <= nextEdgePoint.getRowF(null))
                return PLUS;
            else
                return MINUS;
        }
        else
            return NONE;
    }

    public boolean isRowColAligned()
    {
        return section.isRowColAligned();
    }


    public float[] getRangeZ()
    {
        float zMin;
        float zMax;
        if(points != null)
        {
            zMin = points[0].getZ();
            zMax = zMin;
            for(int n = 1; n < points.length; n++)
            {
                float z = points[n].getZ();
                if(z < zMin) zMin = z;
                else if(z > zMax) zMax = z;
            }
        }
        else
        {
            zMin = prevVertex.getPoint().getZ();
            zMax = zMin;
            float z = nextVertex.getPoint().getZ();
            if(z < zMin) zMin = z;
            else if(z > zMax) zMax = z;
        }
        return new float[] { zMin, zMax };
    }

    public float[] getRangeZorT()
    {
        float zMin;
        float zMax;
        points = getPointsFromEdgePoints();
        if(points != null)
        {
            zMin = points[0].getZorT();
            zMax = zMin;
            for(int n = 1; n < points.length; n++)
            {
                float z = points[n].getZorT();
                if(z < zMin) zMin = z;
                else if(z > zMax) zMax = z;
            }
        }
        else
        {
            zMin = prevVertex.getPoint().getZorT();
            zMax = zMin;
            float z = nextVertex.getPoint().getZorT();
            if(z < zMin) zMin = z;
            else if(z > zMax) zMax = z;
        }
        return new float[] { zMin, zMax };
    }

    static public String getLabel(StsEdge edge)
    {
        if(edge == null) return new String("Null Edge");
        else return edge.getLabel();
    }

    public String labelString()
    {
        String string =  getLabel();
        String prevVertexString = getLabel(prevVertex);
        String nextVertexString = getLabel(nextVertex);
        return string + " from " + prevVertexString + " to " + nextVertexString;
    }

    public void logMessage()
    {
        logMessage(labelString());
    }

    public StsPoint computeEndPointTangent(int end)
    {
        StsPoint point0, point1;

        if(end == MINUS)
        {
            point0 = getFirstPoint();
            point1 = getSecondPoint();
        }
        else if(end == PLUS)
        {
            point0 = getSecondToLastPoint();
            point1 = getLastPoint();
        }
        else
        {
            StsException.outputException(new StsException(StsException.WARNING,
                "StsEdge.getEndPointTangent(end) failed.", " called with: " + end +
                " value must be MINUS(-1) or PLUS(+1)."));
            return null;
        }
        StsPoint point = StsPoint.subPointsStatic(point1, point0);
        return point;
    }

    public void extendEdge(StsPoint point, int direction)
    {
        StsGridSectionPoint gridPoint;
        if(direction == PLUS)
        {
			gridPoint = new StsGridSectionPoint(point);
			nextVertex.setSurfacePoint(gridPoint);
            StsObjectList edgePointsList = getEdgePointsList();
            edgePointsList.add(gridPoint);
        }
        else
        {
			gridPoint = new StsGridSectionPoint(point);
			prevVertex.setSurfacePoint(gridPoint);
			StsList edgePoints = getGridEdgePointsList();
			edgePoints.insertBefore(0, gridPoint);
        }
        points = getPointsFromEdgePoints();
    }

    public void deleteEdgeNotPoints()
    {
        super.delete();
    }

	public void addEdgePoint(StsGridSectionPoint gridPoint)
	{
        StsObjectList edgePointsList = getEdgePointsList();
        if(edgePointsList.getLast() == gridPoint) return;
		edgePointsList.add(gridPoint);
	}

    public StsObjectList getEdgePointsList()
    {
        StsObjectList edgePointsList = null;
        if(isPersistent())
        {
            if(edgePoints == null) edgePoints = StsObjectRefList.constructor(10, 10, "edgePoints", this);
            edgePointsList = edgePoints.getList();
        }
        else
        {
            if(tempEdgePoints == null) tempEdgePoints = new StsObjectList(10, 10);
            edgePointsList = tempEdgePoints;
        }
        return edgePointsList;
    }

    public void addGridEdgePoint(StsGridSectionPoint gridPoint)
	{
		if(gridEdgePoints == null) gridEdgePoints = new StsList(10, 10);
		gridEdgePoints.add(gridPoint);
	}

    public StsPoint findNearestPoint(StsPoint point)
    {
        if(points == null) return null;
        return point.getNearestPointOnSegmentedLine(StsPoint.XY, points, -1.0f);
    }

    public StsPoint findNearestPoint(StsPoint point, float indexF)
    {
        if(points == null) return null;
        return point.getNearestPointOnSegmentedLine(StsPoint.XY, points, indexF);
    }

    public void setPointIndexFs()
    {
        int n;

		if(points == null) return;
        int nPoints = points.length;
        for(n = 0; n < nPoints; n++)
            points[n].setF((float)n);
    }
/*
    protected void insertInterpolatedGridPoint(StsGridSectionPoint gridPoint0, StsGridSectionPoint gridPoint1, float f, StsXYSurfaceGridable grid)
    {
        StsPoint point0, point1, point;

        point0 = gridPoint0.getPoint();
        point1 = gridPoint1.getPoint();
        point = StsPoint.staticInterpolatePoints(point0, point1, f);
        grid.interpolateBilinearZ(point, false, true);
        StsGridSectionPoint gridPoint = new StsGridSectionPoint(point);
        edgePoints.insertAfter(gridPoint0, gridPoint);
        points = getPointsFromEdgePoints();
    }
*/
    StsPoint getTangent(float f)
    {
        int nPoints = points.length;
        int index = (int)((nPoints-1)*f);
        index = StsMath.minMax(index, 0, nPoints-2);
        StsPoint tangent = StsPoint.subPointsStatic(points[index+1], points[index]);
        return tangent;
    }

    public int getConnectedEdgeOrder(int end)
    {
        StsEdge otherEdge = getConnectedEdge(end);
        if(otherEdge == null) return NONE;

        if(end == MINUS)
        {
            if(otherEdge.getPrevVertex() == prevVertex)
                return NOT_OK;
            else
                return OK;
        }
        else  // end == PLUS
        {
            if(otherEdge.getNextVertex() == nextVertex)
                return NOT_OK;
            else
                return OK;
        }
    }

    public StsEdge getConnectedEdge(int end)
    {
        if(end == MINUS)
            return prevVertex.getOtherEdge(this);
        else if(end == PLUS)
            return nextVertex.getOtherEdge(this);
        else
            return null;
    }
/*
    public void reverseConnectedEdges(int end)
    {
        StsEdge edge, nextEdge;

        nextEdge = getConnectedEdge(end);
        while(nextEdge != null)
        {
            edge = nextEdge;
            nextEdge = edge.getConnectedEdge(-end);
            edge.reverseEdge();
        }
    }
*/
    public void reverseEdgePoints()
    {
        if(edgePoints != null)
		{
			edgePoints.reverseOrder();
			points = getPointsFromEdgePoints();
		}
    }

    public void reverseEdge()
    {
        if(edgePoints != null) edgePoints.reverseOrder();
		if(gridEdgePoints != null) gridEdgePoints.reverseOrder();
        if(points != null) StsMath.reverseOrder(points);
		StsSurfaceVertex tempVertex = prevVertex;
		prevVertex = nextVertex;
		nextVertex = tempVertex;
    }

	public StsPoint getPointOnLineNearestMouse(int nSegment, StsMouse mouse, StsGLPanel3d glPanel3d)
	{
		if(points == null) return null;
		int nPoints = points.length;
		int n = Math.min(nPoints-1, nSegment);
        StsMousePoint mousePoint = mouse.getMousePoint();
		return glPanel3d.getPointOnLineNearestMouse(mouse, points[n], points[n+1]);
	}

	// replace the intermediate edgePoints (1 thru nEdgePoints-2) with new edgePoints created from points 1 thru nPoints-2
/*
	public void replaceIntermediateEdgePoints(StsPoint[] points)
	{
		StsGridSectionPoint firstEdgePoint, lastEdgePoint;
		if(edgePoints == null || points == null) return;

		int nEdgePoints = edgePoints.getSize();
		int nPoints = points.length;
		if(nEdgePoints == nPoints) return;

		firstEdgePoint = (StsGridSectionPoint)edgePoints.getFirst();
		lastEdgePoint = (StsGridSectionPoint)edgePoints.getLast();

        for(int n = 1; n < nEdgePoints-1; n++)
        {
            StsGridSectionPoint edgePoint = (StsGridSectionPoint)edgePoints.getElement(n);
            edgePoint.delete();
        }
        edgePoints.delete();
        
        edgePoints = StsObjectRefList.constructor(nPoints, 2, "edgePoints", this);
		edgePoints.add(firstEdgePoint);
		for(int n = 1; n < nPoints-1; n++)
		{
            StsGridSectionPoint edgePoint = new StsGridSectionPoint(points[n], null, true);
            float rowF = section.getRowF(points[n].getZorT());
            float colF = (float)n;
            edgePoint.addSectionRowCol(rowF, colF, section);
//           StsGridSectionPoint edgePoint = new StsGridSectionPoint(points[n], rowF, colF, section, true);
//			edgePoint.setRowOrColIndex(section, COL, n);
			edgePoints.add(edgePoint);
		}
		edgePoints.add(lastEdgePoint);
        constructPoints();
	}
*/
	public void removedDisplayedEdge()
	{
		gridEdgePoints = null;
		points = null;
	}

    public void showPopupMenu(StsWin3d win3d, MouseEvent e)
    {
		if(section == null) return;

        JPopupMenu popup = new JPopupMenu();

        win3d.add(popup);

        StsMenuItem toggleFault = new StsMenuItem();

        toggleFault.setMenuActionListener("Toggle fault", section, "toggleFault", win3d);
		popup.add(toggleFault);

        StsCheckboxMenuItem checkboxSurface = new StsCheckboxMenuItem();

        checkboxSurface.setMenuItemListener("Draw Surface", section, "toggleDrawSurface");
        checkboxSurface.setState(section.isDrawSurface());
		popup.add(checkboxSurface);

        popup.show((Component)win3d, e.getX(), e.getY());
    }
    public String toString()
    {
        String surfaceString = "Surface none";
        if(surface != null)
            surfaceString = surface.toString();
        String sectionString = "Section none";
        if(section != null)
            sectionString = section.toString();
        return sectionString + " " + surfaceString;
    }
}





