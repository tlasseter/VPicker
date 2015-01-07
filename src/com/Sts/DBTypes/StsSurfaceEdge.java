
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

public class StsSurfaceEdge extends StsEdge implements Cloneable, StsEdgeLinkable, StsViewSelectable
{
    protected int side;
	protected int nPointsPerVertex;

    transient private StsBlockGrid blockGrid = null;            /** grid this edge belongs to */
    transient protected StsZoneSide topZoneSide;              /** zoneSide above this edge */
    transient protected StsZoneSide botZoneSide;              /** zoneSide below this edge */
    transient protected float[] arcLengths;                   /** distance along edge: recomputed in several ways */
	transient protected boolean initialized = false;   /** indicates edge added to edgeLoop                */

    static StsSurfaceEdge currentEdge = null;

    // Convenience copies of flags

    static public final int FAULT = StsParameters.FAULT;
    static public final int REFERENCE = StsParameters.REFERENCE;
    static public final int BOUNDARY = StsParameters.BOUNDARY;

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

	public StsSurfaceEdge()
	{
    }

    public StsSurfaceEdge(StsSurfaceVertex v0, StsSurfaceVertex v1,
                          StsModelSurface surface, StsSection section, int side,
                          StsBlockGrid blockGrid)
	{
        super(false);
        this.surface = surface;
        this.section = section;
        this.side = side;
        this.blockGrid = blockGrid;
        addPrevVertex(v0);
        addNextVertex(v1);
//        v0.getSurfacePoint().adjustSectionRowCol(section);
//        v1.getSurfacePoint().adjustSectionRowCol(section);
	}

// Accessors

    public void setSide(int side) { this.side = side; }
    public int getSide() { return side; }
    public StsBlockGrid getBlockGrid() { return blockGrid; }
    public void setTopZoneSide(StsZoneSide topZoneSide) { this.topZoneSide = topZoneSide; }
    public StsZoneSide getTopZoneSide() { return topZoneSide; }
    public void setBotZoneSide(StsZoneSide botZoneSide) { this.botZoneSide = botZoneSide; }
    public StsZoneSide getBotZoneSide() { return botZoneSide; }
    public float[] getArcLengths() { return arcLengths; }

	public StsXYSurfaceGridable getEdgeGrid()
	{
		if(blockGrid != null)
			return blockGrid;
		else
			return surface;
	}

    public StsZoneSide getAZoneSide()
    {
        if(botZoneSide != null) return botZoneSide;
        else return topZoneSide;
    }

    public float getSectionFaultGap()
    {
        if(section == null) return 0.0f;
        return section.getFaultGap(side);
    }

    public StsEdgeLinkable getPrevConnectedSurfaceEdge()
    {
        if(prevVertex == null) return null;
        return prevVertex.getPrevConnectedSurfaceEdge((StsEdgeLinkable)this);
    }


    public StsEdgeLinkable getNextConnectedSurfaceEdge()
    {
        if(nextVertex == null) return null;
        return nextVertex.getNextConnectedSurfaceEdge((StsEdgeLinkable)this);
    }
/*
    public StsGridSectionPoint getOtherSidePoint(int end)
    {
        StsWell well;

        if(end == MINUS)
            well = prevVertex.getWell();
        else
            well = nextVertex.getWell();

        if(well == null) return null;
        return well.getSurfaceEdgeVertex(surface, section, -side, end).getSurfacePoint();
    }
*/
    public boolean onOtherSide(StsBlockGrid otherBlockGrid)
    {
        if(section == null) return false;

        StsList edges = otherBlockGrid.getEdges();
        int nEdges = edges.getSize();

        for(int n = 0; n < nEdges; n++)
        {
            StsSurfaceEdge otherEdge = (StsSurfaceEdge)edges.getElement(n);
            if(otherEdge.getSection() == section && otherEdge.getSide() != side) return true;
        }
        return false;
    }

    public boolean onOtherSide(StsSurfaceEdge otherEdge)
    {
        if(section == null || otherEdge.getSection() != section)
            return false;
        return side != otherEdge.getSide();
    }

    public StsSurfaceEdge overlapsEdgeOnOtherSide(StsBlockGrid otherBlockGrid)
    {
        if(section == null) return null;

        StsList otherEdges = otherBlockGrid.getEdges();
        int nEdges = otherEdges.getSize();

        for(int n = 0; n < nEdges; n++)
        {
            StsSurfaceEdge otherEdge = (StsSurfaceEdge)otherEdges.getElement(n);
            if(otherEdge.getSection() == section && otherEdge.getSide() != side && overlaps(otherEdge)) return otherEdge;
        }
        return null;
    }

    private boolean overlaps(StsSurfaceEdge otherEdge)
    {
        float[] otherColRange = otherEdge.getColRange();
        float[] thisColRange = getColRange();
        return otherColRange[1] > thisColRange[0] && otherColRange[0] < thisColRange[1];
    }

    public boolean isBoundaryOutline()
    {
        if(section == null) return false;
        if(section.getType() != StsParameters.BOUNDARY) return false;
        return side == LEFT;
    }

    public boolean isEndDyingFault(int end)
    {
        if(!isFaulted()) return false;

        StsSurfaceVertex vertex;

        if(end == MINUS)
            return prevVertex.onDyingFault();
        else if(end == PLUS)
            return nextVertex.onDyingFault();
        else
            return false;
    }

    /** draw a segmented dotted-line between points */
    public void display(StsGLPanel3d glPanel3d)
    {
    	StsColor stsColor;
        float edgeWidth;
        boolean debug = false;

		if(section != null)
		{
            int sectionSide = section.getDisplaySide();
            if(side != sectionSide) return;
		}

        StsPoint[] points = getPoints();
        if(points == null) return;

        if(section == null)
        {
            new StsMessage(currentModel.win3d, StsMessage.ERROR, "Section is null\n" +
                           " in StsSurfaceEdge.display()");
            return;
        }

        if(section != null && StsSection.getCurrentSection() == section)
            edgeWidth = StsGraphicParameters.edgeLineWidthHighlighted;
        else
            edgeWidth = StsGraphicParameters.edgeLineWidth;

        GL gl = glPanel3d.getGL();
        if(gl == null) return;

        debug = currentModel.getBooleanProperty("debugDisplayBlockGrids");
        if(!debug)
        {
            if(currentEdge == this)
            {
                stsColor = StsColor.RED;
                prevVertex.display(glPanel3d);
                nextVertex.display(glPanel3d);
            }
            else
                stsColor = StsColor.GREEN;
        }
        else
        {
            if(blockGrid == null)
                stsColor = StsColor.GREY;
            else
            {
                StsSpectrum spectrum = currentModel.getSpectrum("Basic");
                stsColor = spectrum.getColor(blockGrid.getIndex()%32);
            }
        }

        StsGLDraw.drawDottedLine(gl, stsColor, edgeWidth, points);

        if(currentEdge != this) return;

        debug = currentModel.getBooleanProperty("BlockZoneSide Points");
        if(debug)
        {
            if(botZoneSide != null) botZoneSide.getEdgeLoop().debugDisplay(glPanel3d);
            else if(topZoneSide != null) topZoneSide.getEdgeLoop().debugDisplay(glPanel3d);
        }

		if(drawPoints)
		{
			StsGLDraw.drawPoints(points, StsColor.WHITE, glPanel3d,
                                 4, StsGraphicParameters.vertexShift);
		}
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

        if(side == section.getDisplaySide())
            StsGLDraw.pickEdge(gl, color, edgeWidth, points);
    }

    public void mouseSelectedEdit(StsMouse mouse)
    {
        StsBuiltModel builtModel = (StsBuiltModel)currentModel.getCurrentObject(StsBuiltModel.class);
        if(builtModel == null) return;
        builtModel.setPickedSurfaceEdge(this);
        logMessage();
        currentModel.win3dDisplayAll();
    }

    public void showPopupMenu(StsGLPanel glPanel, StsMouse mouse) { }

    /** Make this edge current; set current section. Clear SurfaceEdge if
     *  it is current. Return true if currentEdge was changed.
     */
    public boolean setCurrentEdge()
    {
        currentEdge = this;
        if(section != null) section.setCurrentSection();
        if(blockGrid != null) blockGrid.setCurrentBlockGrid();
        return true;
    }

    static void clearCurrentEdge()
    {
    	currentEdge = null;
    }

    static public StsSurfaceEdge getCurrentEdge()
    {
        return currentEdge;
    }

    /** Should get this directly from Vertex, but we will go thru the
     *  well for now.
     */
    public boolean intersectsVertex(StsSurfaceVertex vertex)
    {
        StsLine line = vertex.getSectionLine();
        if(line == null) return false;

        int lineSectionSide = line.getSectionSide();
        if(lineSectionSide == StsSection.NONE) return false;
        if(lineSectionSide != side) return false;

        StsSection intersectedSection = line.getOnSection();
        if(intersectedSection == section)
        {
            float colF = vertex.getSectionColF(section);
            float colf0 = prevVertex.getSectionColF(section);
            float colf1 = nextVertex.getSectionColF(section);

            return StsMath.betweenInclusive(colF, colf0, colf1);
        }
        return false;
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
        if(blockGrid == null)
            return new String("Edge-" + getIndex() + " on " + section.getLabel());
        else
            return "Edge-" + getIndex() + " on " + section.getLabel() + blockGrid.getLabel();
    }

    public boolean computeArcLengths()
    {
        if(gridEdgePoints == null) return false;
        int nEdgePoints = gridEdgePoints.getSize();
        if(nEdgePoints == 0)
            return false;
        arcLengths = new float[nEdgePoints];

        float length = 0;
        arcLengths[0] = 0.0f;

        StsGridSectionPoint point, nextPoint;

        nextPoint = (StsGridSectionPoint) gridEdgePoints.getFirst();
        for(int n = 1; n < nEdgePoints; n++)
        {
            point = nextPoint;
            nextPoint = (StsGridSectionPoint) gridEdgePoints.getElement(n);
            arcLengths[n] = arcLengths[n-1] + point.getGridDistance(null, nextPoint);
        }
        return true;
    }

    public float getTotalArcLength()
    {
        if(arcLengths == null) return StsParameters.nullValue;
        return arcLengths[arcLengths.length-1];
    }

    // At this point, arc lengths are relative to start of this edge.
    // Adjust so that they are absolute distance to nearest end of this
    // this edge sequence.
    // runningSumLength is the length from the start of the edge sequence up
    // to the start of this edge;
    // totalLength is the length of the total edge sequence
    public float adjustArcLengths(int zeroEnd, float runningSumLength, float totalLength)
    {
        int n;

        if(arcLengths == null) return StsParameters.nullValue;

        int nEdgePoints = arcLengths.length;

        if(runningSumLength > 0.0f)
        {
            for(n = 0; n < nEdgePoints; n++)
                arcLengths[n] += runningSumLength;
        }

        if(zeroEnd == PLUS)
        {
            for(n = 0; n < nEdgePoints; n++)
                arcLengths[n] = totalLength - arcLengths[n];
        }
        else if(zeroEnd == PLUS_AND_MINUS)
        {
            float halfLength = totalLength/2;
            for(n = 0; n < nEdgePoints; n++)
                if(arcLengths[n] > halfLength)
                    arcLengths[n] = totalLength - arcLengths[n];
        }
        return runningSumLength + getTotalArcLength();
    }

    public void setBlockGrid(StsBlockGrid blockGrid)
    {
        this.blockGrid = blockGrid;
    }
}





