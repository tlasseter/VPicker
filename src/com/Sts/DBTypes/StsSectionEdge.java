
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

public class StsSectionEdge extends StsEdge implements Cloneable, StsInstance3dDisplayable, StsViewSelectable
{
	protected int nPointsPerVertex;

	transient static StsSectionEdge currentSectionEdge = null;

    // Convenience copies of flags

    static public final int FAULT = StsParameters.FAULT;
    static public final int REFERENCE = StsParameters.REFERENCE;
    static public final int BOUNDARY = StsParameters.BOUNDARY;
    static public final int AUXILIARY = StsParameters.AUXILIARY;

    static final int MINUS = StsParameters.MINUS;
    static final int PLUS = StsParameters.PLUS;
    static final int NONE = StsParameters.NONE;
    static final int PLUS_AND_MINUS = StsParameters.PLUS_AND_MINUS;

    static final int X = StsPoint.X;
    static final int Y = StsPoint.Y;

    static final float LARGE_FLOAT = StsParameters.largeFloat;
    static final float nullZValue = StsParameters.nullValue;

    public static final int RIGHT = StsParameters.RIGHT;
    public static final int LEFT = StsParameters.LEFT;

    public static final int ROW = StsParameters.ROW;
    public static final int COL = StsParameters.COL;
    public static final int ROWCOL = StsParameters.ROWCOL;

    public static final int FIRST = StsParameters.FIRST;
    public static final int LAST = StsParameters.LAST;

    public static final int OK = StsParameters.OK;
    public static final int NOT_OK = StsParameters.NOT_OK;

// Constructors

	public StsSectionEdge()
	{
    }

	public StsSectionEdge(boolean persistent)
	{
        super(persistent);
    }

    public StsSectionEdge(byte type, StsSection section, StsXYSurfaceGridable surface, int nPointsPerVertex, boolean persistent)
	{
		super(persistent);

        StsLine line;

        this.type = type;
		this.section = section;
		this.surface = surface;
		this.nPointsPerVertex = nPointsPerVertex;

		if(section == null) return;
        StsObjectRefList lines = section.getLines();
        int nLines = lines.getSize();
        edgePoints = StsObjectRefList.constructor(nLines, 2, "edgePoints", this);
        StsGridSectionPoint[] gridSectionPoints = new StsGridSectionPoint[nLines];
        for(int n = 0; n < nLines; n++)
        {
            line = (StsLine)lines.getElement(n);
            StsGridPoint gridPoint = line.computeGridIntersect(surface);
            if(gridPoint == null) continue;
            gridSectionPoints[n] = new StsGridSectionPoint(gridPoint, line);
        }
        edgePoints.add(gridSectionPoints); // add all at once to reduce db write overhead

        line = section.getFirstLine();
        prevVertex = line.getSectionEdgeVertex(surface);
        prevVertex.addEdge(this);
//			convertPointToXYZ(prevVertex);

        line = section.getLastLine();
        nextVertex = line.getSectionEdgeVertex(surface);
        nextVertex.addEdge(this);
    }

	public void convertPointToXYZ(StsSurfaceVertex vertex)
	{
		StsPoint point = vertex.getPoint();
		StsPoint newPoint = point.getXYZorTPoint();
		vertex.setPoint(newPoint);
    }

    public StsSectionEdge(byte type, StsSection section, StsXYSurfaceGridable surface, int nPointsPerVertex)
	{
	    this(type, section, surface, nPointsPerVertex, true);
	}

    public StsSectionEdge(byte type, StsXYSurfaceGridable surface, boolean persistent)
    {
        super(persistent);
        this.type = type;
        this.surface = surface;
    }

    public StsSectionEdge(byte type, StsModelSurface surface)
    {
        this(type, surface, true);
    }

    public boolean initialize(StsModel model)
    {
		StsSeismicVelocityModel velocityModel = currentModel.getProject().getSeismicVelocityModel();
        if (velocityModel != null) adjustTimeOrDepth(velocityModel, section.isOriginalDepth());
        return true;
    }
/*
    public boolean classInitialize()
    {
//		addGridCrossingEdgePoints(true);
		constructPoints();
        return true;
    }
*/
// Accessors

    public void setNPointsPerVertex(int nPointsPerVertex) { this.nPointsPerVertex = nPointsPerVertex; }
    public int getNPointsPerVertex() { return nPointsPerVertex; }

    public boolean isEndDyingFault(int end)
    {
        if(!isFaulted()) return false;

        if(end == MINUS)
            return prevVertex.onDyingFault();
        else if(end == PLUS)
            return nextVertex.onDyingFault();
        else
            return false;
    }
/*
    public void initializeTorZ()
    {
        if(edgePoints == null) return;
        int nEdgePoints = edgePoints.getSize();
        for(int n = 0; n < nEdgePoints; n++)
        {
            StsGridSectionPoint edgePoint = (StsGridSectionPoint)edgePoints.getElement(n);
            StsPoint point = edgePoint.getPoint();
            StsGridPoint gridPoint = new StsGridPoint(point, surface);
            surface.interpolateBilinearZ(gridPoint, true, true);
        }
    }
   */
/*
    public void setPointsFromEdgePoints()
    {
        constructPoints(isDepth);
    }

    public void constructPoints(boolean isDepth)
    {
        StsGridSectionPoint edgePoint;

//		if(points != null) return;

        StsObjectList edgePoints = getEdgePoints();
        if(edgePoints == null) return;

        int nEdgePoints = edgePoints.getSize();
        points = new StsPoint[nEdgePoints];

        for(int i = 0; i < nEdgePoints; i++)
        {
            edgePoint = (StsGridSectionPoint)edgePoints.getElement(i);
            points[i] = edgePoint.getPoint();
//            surface.interpolateBilinearZ(points[i], true, true);
        }

		if(nPointsPerVertex <= 1 || nEdgePoints <= 2) return;
		points = StsBezier.computeXYZorTLine(points, nPointsPerVertex, isDepth);
	}
*/
	public boolean adjustTimeOrDepth(StsSeismicVelocityModel velocityModel, boolean isOriginalDepth)
	{
		StsList edgePoints = getEdgePointsList();
		if(edgePoints == null) return false;
        boolean success = true;
		int nEdgePoints = edgePoints.getSize();
		for(int n = 0; n < nEdgePoints; n++)
		{
			StsGridSectionPoint edgePoint = (StsGridSectionPoint)edgePoints.getElement(n);
			StsPoint point = edgePoint.getPoint();
            if(!velocityModel.adjustTimeOrDepthPoint(point, isOriginalDepth)) success = false;
		}
        return success;
    }

	public boolean isCurved()
	{
		StsList edgePoints = this.getEdgePointsList();
		if(edgePoints == null) return false;
		return edgePoints.getSize() > 2;
	}

    static public void displayClass(StsGLPanel3d glPanel3d, StsClass instanceList)
    {
        boolean displayImportedFaults = true;
//        boolean displayImportedFaults = currentModel.getDisplayImportedFaults();

        int nInstances = instanceList.getSize();
        for(int n = 0; n < nInstances; n++)
        {
            StsSectionEdge edge = (StsSectionEdge)instanceList.getElement(n);
            if(edge.getSection() != null) continue;
            if(displayImportedFaults || edge.getType() != StsParameters.REFERENCE)
                edge.display(glPanel3d);
        }
    }

    public boolean edgeHasSurface(StsXYSurfaceGridable surface)
    {
        return this.surface == surface;
    }
    /** draw a segmented dotted-line between points */
    public void display(StsGLPanel3d glPanel3d)
    {
    	StsColor stsColor;
        float edgeWidth;
		double edgeShift, vertexShift;
        boolean debug = false;

        if(points == null) return;

		// if section built && surface is not visible: don't draw
//        if(!canDrawSection(section)) return;
//        if(section != null && section.getPatch() != null &&
//		   surface != null && !surface.getIsVisible()) return;

        if(section != null && StsSection.getCurrentSection() == section)
            edgeWidth = StsGraphicParameters.edgeLineWidthHighlighted;
        else
            edgeWidth = StsGraphicParameters.edgeLineWidth;

        GL gl = glPanel3d.getGL();
        if(gl == null) return;
/*
        int zIndex;
        if(isDepth)
            zIndex = 2;
        else
            zIndex = 4;
*/
		if(surface != null && surface.getIsVisible())
		{
			edgeShift = StsGraphicParameters.edgeShift;
			vertexShift = StsGraphicParameters.vertexOnEdgeShift;
		}
		else
		{
			edgeShift = 0.0;
			vertexShift = StsGraphicParameters.vertexShift;
		}

        debug = currentModel.getBooleanProperty("debugDisplayBlockGrids");
        if(!debug)
        {
            if(currentSectionEdge == this)
            {
                stsColor = StsColor.RED;
//                prevVertex.display(glPanel3d, StsGraphicParameters.vertexOnEdgeShift, zIndex);
//                nextVertex.display(glPanel3d, StsGraphicParameters.vertexOnEdgeShift, zIndex);
            }
            else if(type == BOUNDARY)
				stsColor = StsColor.CYAN;
			else if(type == FAULT)
                stsColor = StsColor.GREEN;
			else if(type == AUXILIARY)
				stsColor = StsColor.YELLOW;
			else
			    stsColor = StsColor.GREY;
        }
        else
			stsColor = StsColor.GREY;

        gl.glDisable(GL.GL_LIGHTING);

		if(edgeShift != 0.0) glPanel3d.setViewShift(gl, edgeShift);
        StsGLDraw.drawDottedLine(gl, stsColor, edgeWidth, points);
		if(edgeShift != 0.0) glPanel3d.resetViewShift(gl);

		if((drawPoints || currentSectionEdge == this) && points != null && points.length < 20)
		{
			if(nPointsPerVertex == 1)
			{
			    StsGLDraw.drawPoints(points, StsColor.WHITE, glPanel3d, 4, vertexShift);
			}
			else
			{
				int nPoints = points.length;
				for(int n = 0; n < nPoints; n += nPointsPerVertex)
					StsGLDraw.drawPoint(points[n].v, StsColor.WHITE, glPanel3d, 4, vertexShift);
			}
		}
        else if(type != REFERENCE && currentModel.getCurrentObject(StsBuiltModel.class) == null)
        {
            gl.glDisable(GL.GL_LIGHTING);
            checkDisplayVertex(prevVertex, glPanel3d, StsGraphicParameters.vertexOnEdgeShift);
            checkDisplayVertex(nextVertex, glPanel3d, StsGraphicParameters.vertexOnEdgeShift);
            gl.glEnable(GL.GL_LIGHTING);
        }
		gl.glEnable(GL.GL_LIGHTING);
	}

    private boolean canDrawSection(StsSection section)
    {
        if(section == null) return false;
        return section.canDraw();
    }

    private void checkDisplayVertex(StsSurfaceVertex vertex, StsGLPanel3d glPanel3d, double viewShift)
    {
		if(vertex == null) return;
        if(vertex.isFullyConnected()) return;
        vertex.display(glPanel3d, viewShift);
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

//        int zIndex = getPointsZIndex();
        StsGLDraw.pickEdge(gl, color, edgeWidth, points);
    }

    public void mouseSelectedEdit(StsMouse mouse)
    {
        setCurrentEdge();
        if(section != null)
            section.logMessage();
        else
            logMessage();
        currentModel.win3dDisplayAll();
    }

    public void showPopupMenu(StsGLPanel glPanel, StsMouse mouse) { }
/*
    public int getPointsZIndex()
    {
        return currentModel.project.getPointsZIndex();
    }
*/
    /** Make this edge current; set current section. Clear SurfaceEdge if
     *  it is current. Return true if currentEdge was changed.
     */
    public boolean setCurrentEdge()
    {
        currentSectionEdge = this;
        if(section != null) section.setCurrentSection();
        return true;
    }

    static void clearCurrentEdge()
    {
    	currentSectionEdge = null;
    }

    static public StsSectionEdge getCurrentEdge()
    {
        return currentSectionEdge;
    }

    public void setXYZorTPoints(StsPoint[] points)
    {
        int nPoints = points.length;
        this.points = new StsPoint[nPoints];
        for(int n = 0; n < nPoints; n++)
            this.points[n] = points[n].getXYZorTPoint();
    }

    public void showPopupMenu(StsWin3d win3d, MouseEvent e)
    {

        JPopupMenu popup = new JPopupMenu();

        win3d.add(popup);

        StsMenuItem toggleFault = new StsMenuItem();

        toggleFault.setMenuActionListener("Toggle fault", section, "toggleFault", win3d);
		popup.add(toggleFault);

        StsCheckboxMenuItem checkboxSurface = new StsCheckboxMenuItem();

        checkboxSurface.setMenuItemListener("Draw Surface", section, "toggleDrawSurface");
		popup.add(checkboxSurface);
        checkboxSurface.setState(section.isDrawSurface());

        popup.show((Component)win3d, e.getX(), e.getY());
    }

    /** If this edge is on a fault section and this vertex is on a dying fault
     *  (it is not on another fault section): extrapolate edge
     */
    private boolean dyingFaultVertex(StsSurfaceVertex vertex)
    {
        return section.isFault() && vertex.onDyingFault();
    }

    public String getLabel()
    {
		String string = new String("Edge-" + getIndex());
		if(section != null) string = new String(string +  " " + section.getLabel());
		if(surface != null) string = new String(string +  " " + surface.getLabel());
		return string;
    }

    // make a new vertex and split edge at this point; original edge is up to this
    // point and new edge is remaining segment. The two segments each have a vertex
    // at the new point
    public StsSectionEdge[] splitEdge(StsPoint point)
    {
        float vertexF0, vertexF1;
        StsSurfaceVertex vertex;
		StsGridSectionPoint gridPoint;
        StsSectionEdge[] subEdges = new StsSectionEdge[2];

        // this is float index of point location along edge between 0.0 and (float)(nPoints-1)
        // from this, we can compute float index of new vertex location
        // the number of points/vertex is specified by the nPointsPerVertex
        float pointSplitF = point.getF();

		int nEdgePoints = getEdgePointsList().getSize();
        if(pointSplitF < 0.0f || pointSplitF > nEdgePoints-1) return null;

        int prevPointIndex = StsMath.below(pointSplitF);
		gridPoint = new StsGridSectionPoint(point);
        vertex = new StsSurfaceVertex(gridPoint);
        subEdges[0] = makeSubEdge(prevVertex, 0, vertex, prevPointIndex);
        if(subEdges[0] == null) return null;

		gridPoint = new StsGridSectionPoint(point);
        vertex = new StsSurfaceVertex(gridPoint);
        subEdges[1] = makeSubEdge(vertex, prevPointIndex+1, nextVertex, nEdgePoints-1);
        if(subEdges[1] == null) return null;

        deleteEdgeNotPoints();
        return subEdges;
    }

    public StsSectionEdge makeSubEdge(StsSurfaceVertex startVertex, int startIndex, StsSurfaceVertex endVertex, int endIndex)
    {
        StsSectionEdge newEdge = new StsSectionEdge(type, section, surface, nPointsPerVertex);

		if(startVertex != null) // intersection: insert vertex
		    newEdge.addPrevVertex(startVertex);
		else if(prevVertex != null) // no intersection: reuse old vertex
		{
			prevVertex.deleteEdgeFromVertex(this);
			newEdge.addPrevVertex(prevVertex);
		}

		// add intermediate edgePoints
        for(int n = startIndex; n <= endIndex; n++)
        {
            StsGridSectionPoint edgePoint = (StsGridSectionPoint)edgePoints.getElement(n);
            newEdge.addEdgePoint(edgePoint);
        }

		if(endVertex != null) // intersection: insert vertex
		    newEdge.addNextVertex(endVertex);
		else if(nextVertex != null) // no intersection: reuse old vertex
		{
			nextVertex.deleteEdgeFromVertex(this);
			newEdge.addNextVertex(nextVertex);
		}
        newEdge.setPointsFromEdgePoints();
        return newEdge;
    }

    public StsFaultCut convertToFaultCut()
    {
		if(points == null) return null;

        try
        {
		    int nPoints = points.length;
            float[][] localXYs = new float[nPoints][];
            for (int i=0; i<nPoints; i++)
                localXYs[i] = points[i].v;
            return new StsFaultCut(localXYs);
        }
        catch (Exception e)
		{
			StsException.outputException("StsSectionEdge.convertToFaultCut() failed.",
				e, StsException.WARNING);
			return null;
		}
    }

    public void constructFromFaultCut(StsFaultCut faultCut, StsSurface surface)
    {
        if (faultCut == null) return;
        float[][] XYs = faultCut.getXYs();
        if (XYs == null) return;
        int nPnts = XYs.length;

        try
        {
            if(edgePoints == null)
                edgePoints = StsObjectRefList.constructor(nPnts, 2, "edgePoints", this);

            float xMin = (float)(surface.getXMin());
            float xMax = (float)(surface.getXMax());
            float yMin = (float)(surface.getYMin());
            float yMax = (float)(surface.getYMax());

            float minDistSq = 0.25f*surface.getAvgGridIncSq();

            StsPoint point = null, lastPoint;
            for (int i=0; i<nPnts; i++)
            {
                lastPoint = point;
                point = new StsPoint(5);
                float x = StsMath.minMax(XYs[i][0], xMin, xMax);
                point.setX(x);
                float y = StsMath.minMax(XYs[i][1], yMin, yMax);
                point.setY(y);
                surface.interpolateBilinearZ(point, false, true);
                if(i > 0 && i < nPnts-1)
                {
                    float distSq = point.distanceSquared(lastPoint, 3);
                    if(distSq <= minDistSq) continue;
                }

                StsGridSectionPoint edgePoint = new StsGridSectionPoint(point, surface, true);
                edgePoints.add(edgePoint);
            }
            setPointsFromEdgePoints();
            constructEndVertices();
        }
        catch(Exception e)
        {
            StsException.outputException("StsSectionEdge.setVertices() failed.",
                e, StsException.WARNING);
        }
    }

	private void constructEndVertices()
	{
		StsGridSectionPoint edgePoint;

		StsList edgePoints = getEdgePointsList();
		if(edgePoints == null) return;

		edgePoint = (StsGridSectionPoint)edgePoints.getFirst();
	    prevVertex = new StsSurfaceVertex(edgePoint);
		edgePoint = (StsGridSectionPoint)edgePoints.getLast();
		nextVertex = new StsSurfaceVertex(edgePoint);
	}

    // check end vertices: if there is another connected edge check that it is in
    // consistent order - vertex must be the  last vertex for one edge and the first vertex
    // for the next edge.  If one end is out of order and there is not another
    // connected edge at the other end, then reverse this edge; if both ends out of
    // order than reverse this edge; if one end is out of order and there is a connected
    // edge at the other end (which is in order), then all connected edges at the out of
    // order end must be reversed.
/*
    public boolean checkReverseEdgeOrder()
    {
		boolean reversedEdge = false;
        // connection order is OK, NOT_OK (reversed), or NONE (no other connection)
        int minusEndOrder = getConnectedEdgeOrder(MINUS);
        int plusEndOrder = getConnectedEdgeOrder(PLUS);

        if(minusEndOrder == NOT_OK)
        {
            if(plusEndOrder == OK)
                reverseConnectedEdges(MINUS);
            else
			{
                reverseEdge();
				reversedEdge = true;
			}
        }
        else if(plusEndOrder == NOT_OK)
        {
            if(minusEndOrder == OK)
                reverseConnectedEdges(PLUS);
            else
			{
                reverseEdge();
				reversedEdge = true;
			}
        }
		return reversedEdge;
    }
*/
    /*
    public boolean checkEndIntersections(StsXYSurfaceGridable grid)
    {
        if(edgePoints == null || edgePoints.getSize() < 2) return false;
        if(!checkEndIntersection(MINUS, grid)) return false;
        if(!checkEndIntersection(PLUS, grid)) return false;
		return true;
    }
    */
    /*
    private boolean checkEndIntersection(int end, StsXYSurfaceGridable grid)
    {
		StsSurfaceVertex endVertex;
        StsGridSectionPoint endGridPoint, nextGridPoint;
        int startPointIndex, endPointIndex;

        if(points == null) return false;
        int nPoints = points.length;

        if(end == MINUS)
        {
			endVertex = prevVertex;
            endGridPoint = (StsGridSectionPoint)edgePoints.getFirst();
            nextGridPoint = (StsGridSectionPoint)edgePoints.getSecond();
            startPointIndex = 1;
            endPointIndex = nPointsPerVertex-1;
        }
        else
        {
			endVertex = nextVertex;
            endGridPoint = (StsGridSectionPoint)edgePoints.getLast();
            nextGridPoint = (StsGridSectionPoint)edgePoints.getSecondToLast();
            startPointIndex = nPoints - nPointsPerVertex;
            endPointIndex = nPoints-2;
        }

        StsSectionEdge onSectionEdge = endVertex.getSectionEdge();
        if(onSectionEdge == null) return true;
        int correctSide = onSectionEdge.getEdgeSide(nextGridPoint.getPoint());

        // check that intermediate points of segments between first two vertices
        // are on the same side
        for(int n = startPointIndex; n <= endPointIndex; n++)
        {
            int side = onSectionEdge.getEdgeSide(points[n]);
            if(side != correctSide)
            {
                if(end == MINUS)
                    insertInterpolatedGridPoint(endGridPoint, nextGridPoint, 0.99f, grid);
                else
                    insertInterpolatedGridPoint(nextGridPoint, endGridPoint, 0.01f, grid);
                break;
            }
        }
		return true;
    }
    */
    private int getEdgeSide(StsPoint point)
    {
        StsPoint nearestPoint = findNearestPoint(point);
        StsPoint tangent = getTangent(nearestPoint.getF());
        StsPoint normal = StsPoint.subPointsStatic(point, nearestPoint);
        float sideValue = StsPoint.crossProduct2D(normal, tangent);

        if(sideValue > 0.0f) return StsSection.RIGHT;
        else if(sideValue < 0.0f) return StsSection.LEFT;
        else return NONE;
    }

    public void convertFaultToReference(StsWin3d win3d)
    {
        if(type != FAULT) return;
        type = REFERENCE;
        win3d.win3dDisplay();
    }

    public boolean convertFromRefToFault()
    {
        if(type != StsParameters.REFERENCE) return false;
        type = StsParameters.FAULT;
        if(!StsSection.clipEdgeToBoundaries(this)) return false;
        return true;
    }
/*
	public boolean checkAddEndVertices()
	{
		StsGridSectionPoint edgePoint;

		if(edgePoints == null || edgePoints.getSize() < 2) return false;
		if(prevVertex == null)
		{
		    edgePoint = (StsGridSectionPoint)edgePoints.getFirst();
			prevVertex = edgePoint.getVertex();
			if(prevVertex == null)
			    prevVertex = new StsSurfaceVertex((StsGridSectionPoint)edgePoints.getFirst());
			prevVertex.addEdge(this);
		}
		if(nextVertex == null)
		{
		    edgePoint = (StsGridSectionPoint)edgePoints.getLast();
			nextVertex = edgePoint.getVertex();
			if(nextVertex == null)
			    nextVertex = new StsSurfaceVertex((StsGridSectionPoint)edgePoints.getLast());
			nextVertex.addEdge(this);
		}
		return true;
	}
*/
	public void addEdgeToConnectedVertices()
	{
		StsSurfaceVertex[] connectedVertices = getConnectedVertices();
		if(connectedVertices == null) return;
		int nConnectedVertices = connectedVertices.length;
		for(int n = 0; n < nConnectedVertices; n++)
			connectedVertices[n].addEdgeAssociation(this);
	}

	private StsSurfaceVertex[] getConnectedVertices()
	{
		StsSurfaceVertex[] connectedVertices = null;

		StsList edgePoints = getEdgePointsList();
		if(edgePoints == null) return null;
		int nEdgePoints = edgePoints.getSize();
		if(nEdgePoints <= 2) return null;
		for(int n = 1; n < nEdgePoints-1; n++)
		{
			StsGridSectionPoint edgePoint = (StsGridSectionPoint)edgePoints.getElement(n);
			StsSurfaceVertex vertex = edgePoint.getVertex();
			if(vertex == null) continue;
			connectedVertices = (StsSurfaceVertex[])StsMath.arrayAddElement(connectedVertices, vertex);
		}
		return connectedVertices;
	}

	public boolean delete()
	{
		removeEdgeFromConnectedVertices();
		super.delete();
        currentModel.removeDisplayableInstance(this);
        return true;
	}

	private void removeEdgeFromConnectedVertices()
	{
		StsSurfaceVertex[] connectedVertices = getConnectedVertices();
		if(connectedVertices == null) return;
		int nConnectedVertices = connectedVertices.length;
		for(int n = 0; n < nConnectedVertices; n++)
			connectedVertices[n].deleteEdge(this);
	}

	public StsSection[] getConnectedSections()
	{
		return null;
	}

	public StsLine[] getLines()
	{
		StsLine[] lines = new StsLine[2];
		if(prevVertex != null && nextVertex != null)
		{
		    lines[0] = prevVertex.getSectionLine();
			lines[1] = nextVertex.getSectionLine();
		}
		return lines;
	}

}





