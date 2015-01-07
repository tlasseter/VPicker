
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
import com.Sts.Utilities.Interpolation.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.util.*;

/** StsBlockGrid is a persistent grid structure defining the trimmed-grid
 *  resulting from the intersection of a horizon with a block
 */

public class StsBlockGrid extends StsGridBoundingBox implements Cloneable, StsXYSurfaceGridable
{

    // instance fields
    transient protected StsModelSurface surface;
    transient public StsBlock block;
    transient protected StsEdgeLoop edgeLoop;
    transient protected StsBlockGrid blockGridBelow;
    transient protected StsBlockGrid blockGridAbove;

    transient protected int pointRowMin, pointRowMax, pointColMin, pointColMax;  // Non-null entry limits
	transient protected int nRows, nCols;
    transient protected float angle = 0.0f;
    transient protected float nullZValue;

    transient protected float[][] pointsZ = null;
    transient protected byte[][] pointsType = null;

    transient protected float[][] thickness = null;
    transient protected float[][][] normals = null;

    transient int surfaceNRows, surfaceNCols;
    transient protected int newPointRowMin, newPointRowMax, newPointColMin, newPointColMax; // temps used in resetting grid point limits

	transient StsObjectList connectedBlockGrids; // block grid adjoing this one across auxiliary section

    transient StsCellTypeGrid cellTypeGrid;  // used only for debugging

    static StsBlockGrid currentBlockGrid = null;
    static protected StsBlockGridInterpolation interpolator;

    static boolean doDebugInterpolate = true; // we want to mainDebug interpolation during construction
    static boolean debugInterpolate;
    static int debugInterpolateRow = 56;
    static int debugInterpolateCol = 1;

    static byte debugInterpolatePointType;

//    transient StsList quadStrips;
//    transient StsList polygons;

//    transient boolean gridConstructed = false;

//    transient StsLinkGridLine[] rowGridLines;
//    transient StsLinkGridLine[] colGridLines;

    static final float minDeterminant = 1.0e-5f;

//    transient boolean useGapOrGrid = false; // if false, interpolate gapPoint if null
                                            // if true, use surface point if gapPoint is null

    transient float[] vertNormal = new float[] { 0.0f, 0.0f, -1.0f };
    static float[] drawPoint = new float[3];

//    transient StsWin3d win3d;  // for mouse & message handling

    // pointsType flags
    static public final byte GAP_NULL = StsParameters.GAP_NULL;             // Not initialized
    static public final byte GAP_GRID = StsParameters.GAP_GRID;             // gridPoint for interpolating gapPoints
    static public final byte GAP_SURF_GRID = StsParameters.GAP_SURF_GRID;   // surface gridPoint (used only in interpolating unfaulted edges) not used
    static public final byte GAP_NOT_FILLED = StsParameters.GAP_NOT_FILLED; // gapPoint we need to fill
    static public final byte GAP_FILLED = StsParameters.GAP_FILLED;         // filled by least-sqs interpolated value
    static public final byte GAP_CANT_FILL = StsParameters.GAP_CANT_FILL;   // couldn't compute least-sq value
    static public final byte GAP_CUT = StsParameters.GAP_CUT;               // gapPoint, but has been cutoff by fault
    static public final byte GAP_EXTRAP = StsParameters.GAP_EXTRAP;         // extrapolated by least-sqs interpolated value
    static public final byte GAP_NONE = StsParameters.GAP_NONE;             // passed when type not to be changed
    static public final byte GAP_FILL_HOLE = StsParameters.GAP_FILL_HOLE;   // back-filled hole between edge and first GAP_GRID point: not used

    // Convenience copies of useful flags
    static public final int PLUS = StsParameters.PLUS;
    static public final int MINUS = StsParameters.MINUS;
    static public final float nullValue = StsParameters.nullValue;
    static public final float roundOff = StsParameters.roundOff;
    static public final float largeFloat = StsParameters.largeFloat;

    static public final int ROW = StsParameters.ROW;
    static public final int COL = StsParameters.COL;
    static public final int ROWCOL = StsParameters.ROWCOL;

    static public final int NONE = StsParameters.NONE;

    static public final byte QUAD = 1;
    static public final byte POLY = 2;
    static public final byte EMPTY = 3;

    static public final int ROW_PLUS = StsParameters.ROW_PLUS;
    static public final int COL_PLUS = StsParameters.COL_PLUS;
    static public final int ROW_MINUS = StsParameters.ROW_MINUS;
    static public final int COL_MINUS = StsParameters.COL_MINUS;

    static public final int nullInteger = StsParameters.nullInteger;

    /** Convenience copies of TStrip flags */
    static public final byte STRIP_INVALID = StsParameters.STRIP_INVALID;
    static public final byte STRIP_BOTH = StsParameters.STRIP_BOTH;
    static public final byte STRIP_BOT = StsParameters.STRIP_BOT;
    static public final byte STRIP_TOP = StsParameters.STRIP_TOP;

    static public final byte GAP_LINE = StsParameters.GAP_LINE;
    static public final byte GRID_LINE = StsParameters.GRID_LINE;
    static public final byte BOUNDARY_LINE = StsParameters.BOUNDARY_LINE;

    static public final int GRIDLINE_GRID = StsParameters.GRIDLINE_GRID;
    static public final int GRIDLINE_BOUNDARY = StsParameters.GRIDLINE_BOUNDARY;
    static public final int GRIDLINE_NONE = StsParameters.GRIDLINE_NONE;

    static public final byte FINISHED = 1;
    static public final byte FILLED = 2;
    static public final byte CONTINUE = 3;

    static final byte CELL_EMPTY = StsParameters.CELL_EMPTY;
    static final byte CELL_FILL = StsParameters.CELL_FULL;
    static final byte CELL_EDGE = StsParameters.CELL_EDGE;

    static final boolean SET_POINT_TRUE = true;
    static final boolean SET_POINT_FALSE = false;
    static final boolean DEBUG_TRUE = true;
    static final boolean DEBUG_FALSE = false;
    static final boolean INSIDE_TRUE = true;
    static final boolean INSIDE_FALSE = false;

	public StsBlockGrid()
	{
	}

    public StsBlockGrid(StsModelSurface surface, StsBlock block)
    {
        super(false);
        this.surface = surface;
        this.block = block;

        if(block != null)
        {
            edgeLoop = new StsEdgeLoop(this);
            edgeLoop.setNormal(vertNormal, false);
        }

        surfaceNRows = surface.getNRows();
        surfaceNCols = surface.getNCols();

        newPointRowMin = surfaceNRows-1;
        newPointRowMax = 0;
        newPointColMin = surfaceNCols-1;
        newPointColMax = 0;
    }

    public boolean constructInitialEdges()
    {
        StsObjectRefList blockSides;
        StsBlockSide blockSide;
        StsSurfaceEdge surfaceEdge;

        try
        {
            blockSides = block.getBlockSides();
            int nBlockSides = blockSides.getSize();

            for(int n = 0; n < nBlockSides; n++)
            {
                blockSide = (StsBlockSide)blockSides.getElement(n);
                surfaceEdge = blockSide.constructInitialEdge(surface, this);
                edgeLoop.addEdge(surfaceEdge, PLUS);
            }
            constructGridBoundingBox();
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsBlockGrid.constructEdges() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    /** construct an XYZ loopBoundingBox and set the row and col ranges. */
    public void constructGridBoundingBox()
    {
        edgeLoop.computeGridBoundingBoxFromEdges(this);
        integerAdjustRange(surface);
        computeRowColRanges(surface);
    }

    public void initializeArrays()
    {
        try
        {
            pointRowMin = block.rowMin;
            pointRowMax = block.rowMax;
            pointColMin = block.colMin;
            pointColMax = block.colMax;
            addPointBorder(2, surface);
/*
            int border = 1;

            pointRowMin = Math.max(newPointRowMin-border, 0);
            pointRowMax = Math.min(newPointRowMax+border, surfaceNRows-1);
            pointColMin = Math.max(newPointColMin-border, 0);
            pointColMax = Math.min(newPointColMax+border, surfaceNCols-1);
*/
            nRows = pointRowMax - pointRowMin + 1;
            nCols = pointColMax - pointColMin + 1;

            if(nRows <= 0 || nCols <= 0)
            {
                nRows= 2; nCols = 2;
                pointRowMin = 0; pointRowMax = 1;
                pointColMin = 0; pointColMax = 1;
            }

            newPointRowMin = surfaceNRows-1;
            newPointRowMax = 0;
            newPointColMin = surfaceNCols-1;
            newPointColMax = 0;

            xInc = surface.xInc;
            yInc = surface.yInc;

            // this is only approximate, but we don't need accurate limits (yet)
            zMin = surface.getZMin();
            zMax = surface.getZMax();

            nullZValue = surface.nullZValue;

            initPointsZ();
            initPointsType();
        }
        catch(Exception e)
        {
            StsException.outputException("StsBlockGrid.initializeArrays() failed.",
                e, StsException.WARNING);
        }
    }

    public boolean initialize(StsModel model)
    {
        return construct();
    }

    public void addToLinkedGrid(StsEdgeLoopRadialLinkGrid linkedGrid)
    {
        if(!edgeLoop.constructEdgeLinks()) return;
        // edgeLoop.addToLinkedGrid(linkedGrid);
        linkedGrid.addEdgeLoop(edgeLoop);
    }

    public void construct2dGrid(StsModel model)
    {
		setConnectedBlockGrids();
        checkInterpolator();
        fillInsideGridPoints();
        intersectGridWithEdges(model);
        adjustGridRangeForAuxiliaryEdges();
    }

    private void checkInterpolator()
    {
        if(interpolator == null)
            interpolator = StsBlockGridInterpolationWeightedPlane.getInstance(surface, this);
        interpolator.initializeGrid(surface, this);
//        interpolator = StsBlockGridInterpolationLstSqPlaneFit.getInstance(surface, this);
    }

    private void setConnectedBlockGrids()
	{
		StsObjectRefList connectedBlocks = block.getConnectedBlocks();
		if(connectedBlocks == null) return;
		int nConnectedBlocks = connectedBlocks.getSize();
		if(nConnectedBlocks == 0) return;
		connectedBlockGrids = new StsObjectList(nConnectedBlocks);
		for(int n = 0; n < nConnectedBlocks; n++)
		{
			StsBlock block = (StsBlock)connectedBlocks.getElement(n);
			StsBlockGrid connectedBlockGrid = block.getBlockGrid(surface);
			if(connectedBlockGrid != null) connectedBlockGrids.add(connectedBlockGrid);
		}
	}

	public boolean connectsBlockGrid(StsBlockGrid otherBlockGrid)
	{
		if(connectedBlockGrids == null) return false;
		int nConnectedBlockGrids = connectedBlockGrids.getSize();
		for(int n = 0; n < nConnectedBlockGrids; n++)
		{
			if(otherBlockGrid == connectedBlockGrids.getElement(n))
				return true;
		}
		return false;
	}

    public boolean construct()
    {
        if(edgeLoop == null) return false;
        return edgeLoop.construct(this);
    }

    private void intersectGridWithEdges(StsModel model)
    {
        StsSurfaceEdge surfaceEdge;
        int n;

        StsList edges = getEdges();
        int nEdges = edges.getSize();

        for(n = 0; n < nEdges; n++)
        {
            surfaceEdge = (StsSurfaceEdge)edges.getElement(n);
            surfaceEdge.getNextVertex().adjustToGrid(this);
        }

        for(n = 0; n < nEdges; n++)
        {
            surfaceEdge = (StsSurfaceEdge)edges.getElement(n);
            surfaceEdge.constructSurfaceEdgePoints(this, false);
        }
    }

    // Accessors
    public float getAngle() { return angle; }
    public int getNRows() { return nRows; }
    public int getNCols() { return nCols; }
//    public float getXInc() { return xInc; }
//    public float getYInc() { return yInc; }
    public float getZMin() { return zMin; }
    public float getZMax() { return zMax; }
    public int getPointRowMin() { return pointRowMin; }
    public int getPointRowMax() { return pointRowMax; }
    public int getPointColMin() { return pointColMin; }
    public int getPointColMax() { return pointColMax; }
    public int getRowMin() { return pointRowMin; }
    public int getRowMax() { return pointRowMax; }
    public int getColMin() { return pointColMin; }
    public int getColMax() { return pointColMax; }
    public StsModelSurface getSurface() { return surface; }
    public StsEdgeLoop getEdgeLoop() { return edgeLoop; }
    public StsGridBoundingBox getGridBoundingBox() { return this; }
    public boolean hasAuxiliaryEdges() { return block.hasAuxiliarySections(); }
    public void setBlockGridAbove(StsBlockGrid blockGrid) { blockGridAbove = blockGrid; }
    public void setBlockGridBelow(StsBlockGrid blockGrid) { blockGridBelow = blockGrid; }
    public StsBlockGrid getBlockGridAbove() { return blockGridAbove; }
    public StsBlockGrid getBlockGridBelow() { return blockGridBelow; }
    public StsColor getStsColor() { return surface.getStsColor(); }

    public double getXOrigin() { return surface.getXOrigin(); }
    public double getYOrigin() { return surface.getYOrigin(); }
    public float getXSize() { return surface.getXSize(); }
    public float getYSize() { return surface.getYSize(); }
    public float getXMin(){ return surface.getXMin(); }
    public float getYMin() { return surface.getYMin(); }
    public float getXInc() { return surface.getXInc(); }
    public float getYInc() { return surface.getYInc(); }
    public float getZInc() { return 0.0f; } // for StsXYSurfaceGridable compatability
    public float getRowCoor(float[] xy) { return surface.getRowCoor(xy); }
    public float getColCoor(float[] xy) { return surface.getColCoor(xy); }
    public float getYCoor(float rowF, float colF) { return surface.getYCoor(rowF, colF); }
    public float getXCoor(float rowF, float colF) { return surface.getXCoor(rowF, colF); }


    public StsSurfaceEdge getSurfaceEdge(StsSection section)
    {
        StsList edgeList = getEdges();
        int nEdges = edgeList.getSize();
        for(int n = 0; n < nEdges; n++)
        {
            StsEdgeLinkable edge = (StsEdgeLinkable)edgeList.getElement(n);
            if(!(edge instanceof StsSurfaceEdge)) continue;
            StsSurfaceEdge surfaceEdge = (StsSurfaceEdge)edge;
            if(surfaceEdge.section == section) return surfaceEdge;
        }
        return null;
    }

    public StsList getEdges()
    {
        if(edgeLoop == null) return new StsList(0);
        else return edgeLoop.getEdges();
    }

    public void setPointType(int row, int col, byte type)
    {
        try
        {
            byte currentType = pointsType[row-pointRowMin][col-pointColMin];

            if(type == currentType) return;

            adjustRowColRange(row, col);

            if(Main.debugPoint && row == Main.debugPointRow && col == Main.debugPointCol)
            {
                StsException.systemDebug("StsBlockGrid.setPointType() called. " +
                    this.getLabel() + "pointsType[" + row + "][" + col + "] changed to " +
                    StsParameters.getGapTypeName(type) + " from " +
                    StsParameters.getGapTypeName(currentType) );
            }

            /** Type changes are complicated, but can best be understood in sequence. */

            switch(type)
            {
                /** Under some circumstances, we reset the pointType to NULL so that it
                 *  can be subsequently recomputed and reset to another type.
                 */
                case GAP_NULL:
                {
                    pointsType[row-pointRowMin][col-pointColMin] = type;
                    break;
                }
                /** Initially, pointType is NULL; points are set to GRID meaning they are
                 *  good points that can be used for interpolation or GAP_NOT_FILLED meaning
                 *  they are in a fault gap and need to be filled.
                 */
                case GAP_GRID:
                case GAP_SURF_GRID:
                case GAP_NOT_FILLED:
                case GAP_CANT_FILL:
                {
                    pointsType[row-pointRowMin][col-pointColMin] = type;
                    setPointZ(row, col, surface.getZorT(row, col), type);
                    break;
                }

                /** Least-squares interpolation/extrapolation is used to fill points need for
                 *  blockGrid.  GAP_NOT_FILLED points are filled and flagged as GAP_FILLED or
                 *  GAP_CANT_FILL (happens only if there are too few points in which case they
                 *  retain the original grid point.  When the edge is intersected against this
                 *  blockGrid, we may need other extrapolated points which are filled and flagged
                 *  with GAP_EXTRAP.
                 */
                case GAP_FILLED:
                case GAP_EXTRAP:
                case GAP_FILL_HOLE:
                {
                    pointsType[row-pointRowMin][col-pointColMin] = type;
                    break;
                }

                /** After intersected edges have been computed, we can remove any points
                 *  which are deemed to be on the wrong side of the intersection. Flag
                 *  these points as GAP_CUT (CAP_GUT?)
                 */

                case GAP_CUT:
                {
                    pointsType[row-pointRowMin][col-pointColMin] = type;
                    break;
                }

                /** Passed to leastSqFitInterpolate indicating we want the value but don't
                 *  change the type.  Used in determining surface position when a value needs
                 *  to be extrapolated in order to compute, but we don't want the type changed.
                 */
                case GAP_NONE:
                {
                    break;
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsBlockGrid.setPointType() failed.",
                e, StsException.WARNING);
        }
    }

    public void setPointTypeIfNot(int row, int col, byte type, byte notType)
    {
        if(isPointType(row, col, notType)) return;
        setPointType(row, col, type);
    }

    public void setPointTypeIf(int row, int col, byte type, byte ifType)
    {
        if(!isPointType(row, col, ifType)) return;
        setPointType(row, col, type);
    }

    // we need to include the range of the auxiliary edge points because the grid
    // we subsequently be downsized and these points need to be included
    public void adjustGridRangeForAuxiliaryEdges()
    {
        StsEdgeLoop edgeLoop = getEdgeLoop();
        if(edgeLoop == null) return;

        StsList edges = getEdges();
        int nEdges = edges.getSize();
        for(int e = 0; e < nEdges; e++)
        {
            StsSurfaceEdge edge = (StsSurfaceEdge)edges.getElement(e);
            if(!edge.isAuxiliary()) continue;

            Iterator<StsEdgeLoopRadialGridLink> edgeLinksIterator = edgeLoop.getEdgeLinksIterator(edge, PLUS);
            while(edgeLinksIterator.hasNext())
            {
                StsEdgeLoopRadialGridLink link = edgeLinksIterator.next();
                int[] rowCol = link.getOutsideAdjacentRowAndCol();
                if(rowCol != null) adjustRowColRange(rowCol[0], rowCol[1]);
            }
        }
    }

    public void adjustRowColRange(int row, int col)
    {
        newPointRowMin = Math.min(newPointRowMin, row);
        newPointRowMax = Math.max(newPointRowMax, row);
        newPointColMin = Math.min(newPointColMin, col);
        newPointColMax = Math.max(newPointColMax, col);
    }

    public byte getPointType(int row, int col)
    {
        if(pointsType == null) return GAP_NULL;
        if(insideGrid(row, col)) return pointsType[row-pointRowMin][col-pointColMin];
//        StsException.systemError("StsBlockGrid.getPointType() failed." + outsideGridMessage(row, col));
        return GAP_NULL;
    }

    public void setPointsType(int rowOrCol, int index, int min, int max, byte type)
    {
        if(rowOrCol == ROW)
            for(int j = min; j <= max; j++)
                setPointType(index, j, type);
        else // rowOrCol == COL
            for(int i = min; i <= max; i++)
                setPointType(i, index, type);
    }

    public boolean isPointType(int row, int col, byte type)
    {
        return getPointType(row, col) == type;
    }

    public boolean changePointTypeIfNot(int row, int col, byte excludeType, byte newType)
    {
        if(getPointType(row, col) == excludeType) return false;
        setPointType(row, col, newType);
        return true;
    }

    public boolean changePointTypeIf(int row, int col, byte changeType, byte newType)
    {
        if(getPointType(row, col) != changeType) return false;
        setPointType(row, col, newType);
        return true;
    }

    /** Called by StsLineLinkGrid to determine if this is a gap point */
    public boolean isPointGap(int i, int j)
    {
        if(pointsType == null) return false;
        return isPointGap(getPointType(i, j));
    }

    public boolean isPointGapDomain(int i, int j)
    {
        if(pointsZ == null) return false;
        return isPointGapDomain(getPointType(i, j));
    }
/*
    public boolean isPointInDomain(int i, int j)
    {
        if(pointsZ == null) return false;
        if(!insideCellGrid(i, j)) return false;
        return isPointInDomain(getPointType(i, j));
    }
*/
    public boolean isPointNotActive(int i, int j)
    {
        return isPointNotActive(getPointType(i, j));
    }

    public boolean isPointNotActive(byte pointType)
    {
        return pointType == GAP_NULL ||
               pointType == GAP_CANT_FILL ||
               pointType == GAP_NOT_FILLED ||
               pointType == GAP_CUT;
    }

    static public boolean isPointGap(byte pointType)
    {
        return pointType != GAP_GRID && pointType != GAP_SURF_GRID;
    }

    final public boolean isPointGrid(int i, int j)
    {
        if(pointsZ == null) return false;
        return isPointGrid(getPointType(i, j));
    }

    static public boolean isPointGrid(byte pointType)
    {
        return pointType == GAP_GRID;
//        return pointType == GAP_GRID || pointType == GAP_SURF_GRID;
    }

    public boolean isPointGridOrSurface(int i, int j)
    {
        if(pointsZ == null) return false;
        return isPointGridOrSurface(getPointType(i, j));
    }

    static public boolean isPointGridOrSurface(byte pointType)
    {
        return pointType == GAP_GRID || pointType == GAP_SURF_GRID;
    }

    static public boolean isPointGapDomain(byte pointType)
    {
        return pointType == GAP_FILLED ||
               pointType == GAP_CANT_FILL ||
               pointType == GAP_FILL_HOLE;
    }

    static public boolean isPointInDomain(byte pointType)
    {
        return pointType == GAP_GRID ||
               pointType == GAP_SURF_GRID ||
               isPointGapDomain(pointType);
    }

    static public boolean isFillable(byte pointType)
    {
        return pointType == GAP_NULL || pointType == GAP_NOT_FILLED || pointType == GAP_SURF_GRID || pointType == GAP_EXTRAP;
    }

    // StsBlockGrid should go away; set & getProperty is here for interface compatability
    public void setProperty(int row, int col, float value) { return; }
    public float getProperty(int row, int col) { return StsParameters.nullValue; }

    public void setCurrentBlockGrid() { currentBlockGrid = this; }
    static public StsBlockGrid getCurrentBlockGrid() { return currentBlockGrid; }

    public float[][] getPointsZ() { return pointsZ; }

    public byte[][] getPointsType() { return pointsType; }

    public final float[] getXYZorT(int row, int col)
    {
        float x, y, z;
        z = getPointZ(row, col);
        if(z == nullZValue) return null;
        x = surface.getXCoor(row, col);
        y = surface.getYCoor(row, col);
        return new float[] { x, y, z };
    }

    public StsPoint getPoint(int row, int col)
    {
        float x, y, z;
        z = getPointZ(row, col);
        if(z == nullZValue) return null;
        x = surface.getXCoor(row, col);
        y = surface.getYCoor(row, col);
        return new StsPoint( x, y, z );
    }

    public StsGridPoint getGridPoint(int row, int col)
    {
        if(pointsZ == null) return null;
        if(getPointType(row, col) != GAP_GRID) return null;

        StsGridPoint gridPoint = new StsGridPoint(this);
        gridPoint.setRowCol(row, col);
        float z = getPointZ(row, col);
        gridPoint.setZ(z);
        gridPoint.setT(z);
        gridPoint.setNullType(StsParameters.SURF_GAP);
        return gridPoint;
    }

    public StsGridPoint getGapPoint(int row, int col)
    {
        if(pointsZ == null || !isPointGap(row, col)) return null;

        StsGridPoint gridPoint = new StsGridPoint(this);
        gridPoint.setRowCol(row, col);
        float z = getPointZ(row, col);
        gridPoint.setZ(z);
        gridPoint.setT(z);
        gridPoint.setNullType(StsParameters.SURF_GAP);
        return gridPoint;
    }

    public StsGridPoint getGridOrGapPoint(int row, int col)
    {
        if(pointsZ == null) return null;

        byte pointType = getPointType(row, col);

        if(isPointGap(row, col))
        {
            StsGridPoint gridPoint = new StsGridPoint(this);
            gridPoint.setRowCol(row, col);
            float z = getPointZ(row, col);
            gridPoint.setZ(z);
            gridPoint.setT(z);
            gridPoint.setNullType(StsParameters.SURF_GAP);
            return gridPoint;
        }
        else
            return surface.getGridPoint(row, col);
    }

    public float getGridOrGapPointZ(int row, int col)
    {
        if(pointsZ == null) return nullZValue;

        byte pointType = getPointType(row, col);

        if(isPointGap(row, col))
            return getPointZ(row, col);
        else
            return surface.getZorT(row, col);
    }

    public float[] getXYZorT(float rowF, float colF)
    {
        StsPoint point = getPoint(rowF, colF);
        return point.getXYZorT();
    }
    public StsPoint getPoint(float rowF, float colF)
    {
        try
        {
            int row = (int)rowF;
            float dR = rowF - row;

            int col = (int)colF;
            float dC = colF - col;

            if(dR == 0.0f)
            {
                if(dC == 0.0f)
                    return getPoint(row, col);
                else
                    return StsPoint.staticInterpolatePoints(getPoint(row, col), getPoint(row, col+1), dC);
            }
            else
            {
                StsPoint point0 = StsPoint.staticInterpolatePoints(getPoint(row, col), getPoint(row+1, col), dR);
                if(dC == 0.0f)
                    return point0;
                else
                {
                    StsPoint point1 = StsPoint.staticInterpolatePoints(getPoint(row, col+1), getPoint(row+1, col+1), dR);
                    return StsPoint.staticInterpolatePoints(point0, point1, dC);
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsBlockGrid.getNormal() failed.", e, StsException.WARNING);
            return null;
        }
    }

   /* check for z null */
    public boolean okGridPoint(int row, int col)
    {
        return getPointType(row, col) == GAP_GRID;
    }

    private void setNormalXYZ(int row, int col, float[] normal)
    {
        if(getPointType(row, col) == GAP_NULL) return;
        normals[row][col] = new float[] { normal[0], normal[1], normal[2] };
    }

    public boolean insideGrid(int row, int col)
    {
        return row >= pointRowMin && row <= pointRowMax && col >= pointColMin && col <= pointColMax;
    }

    public boolean insideGrid(StsGridPoint gridPoint)
    {
        return insideGrid(gridPoint.row, gridPoint.col);
    }

    public String outsideGridMessage(int row, int col)
    {
        if(row < pointRowMin) return new String(" row: " + row + " < " + pointRowMin);
        if(row > pointRowMax) return new String(" row: " + row + " > " + pointRowMax);
        if(col < pointColMin) return new String(" col: " + col + " < " + pointColMin);
        if(col > pointColMax) return new String(" col: " + col + " > " + pointColMax);
        return new String(" ");
    }

    public boolean insideGrid(float rowF, float colF)
    {
        return rowF >= pointRowMin && rowF <= pointRowMax && colF >= pointColMin && colF <= pointColMax;
    }

    public void addPointBorder(int border, StsXYGridable limitGrid)
    {
        pointRowMin = Math.max(limitGrid.getRowMin(), pointRowMin - border);
        pointRowMax = Math.min(limitGrid.getRowMax(), pointRowMax + border);
        pointColMin = Math.max(limitGrid.getColMin(), pointColMin - border);
        pointColMax = Math.min(limitGrid.getColMax(), pointColMax + border);
    }
    
    private void initPointsZ()
    {
        pointsZ = new float[nRows][nCols];

        for(int i = 0; i < nRows; i++)
            for(int j = 0; j < nCols; j++)
                pointsZ[i][j] = nullZValue;
    }

    private void initPointsType()
    {
       pointsType = new byte[nRows][nCols];
/*
        for(int i = pointRowMin; i <= pointRowMax; i++)
        {
            for(int j = pointColMin; j <= pointColMax; j++)
            {
                {
                    switch(surface.getPointNull(i, j))
                    {
                        case StsParameters.NOT_NULL:
                            setPointType(i, j, GAP_SURF_GRID);
                            break;
                        case StsParameters.NULL_GAP:
                        case StsParameters.NULL_GAP_SET:
                    }
                }
            }
        }
*/
    }
/*
    public void intersectAgainstGridBelow()
    {
        if(blockGridBelow == null) return;

        for(int i = pointRowMin; i <= pointRowMax; i++)
            for(int j = pointColMin; j <= pointColMax; j++)
                if(!isPointGap(i, j)) intersectAgainstGridBelow(i, j);
    }

    public void intersectAgainstGridBelow(int row, int col)
    {
        float z = getPointZ(row, col);
        setPointZ(row, col, z, fillFlag);
    }
*/
    public void constructLeastSqBlockGrid()
    {
        for(int i = pointRowMin; i <= pointRowMax; i++)
            for(int j = pointColMin; j <= pointColMax; j++)
                if(getPointType(i, j) == GAP_NOT_FILLED)
                    interpolateInsidePoint(i, j);
    }

    public void constructLeastSqBlockGridWithThickness(StsModelSurface baseSurface, StsModelSurface thickness)
    {
    /*
        for(int i = 0; i < nRows; i++)
        {
            for(int j = 0; j < nCols; j++)
            {
                float z = leastSqFitGridPoint(i, j, baseSurface);
                if(z != nullZValue) pointsZ[i][j] = z - thickness.getPointZ(i, j);
            }
        }
    */
    }

    // If we know a point is inside the edgeLoop and there are auxiliary edges,
    // then we can weight points outside the loop
    public float interpolateInsidePoint(int in, int jn)
    {
        return interpolate(in, jn, INSIDE_TRUE, SET_POINT_TRUE, DEBUG_FALSE);
    }

    // If this point is outside the loop, we do not want to use any other points
    // outside the loop for interpolation.
    // If setPoint is false, we do not set the point value and type.
    public float interpolateOutsidePoint(int in, int jn, boolean setPoint)
    {
        return interpolate(in, jn, INSIDE_FALSE, setPoint, DEBUG_FALSE);
    }

    // Spiral around this point, looking for good points to interpolate.  A point is good if it is
    // inside the edgeLoop and of type GAP_GRID or this point itself is inside (isInside == true),
    // the edgeLoop has auxiliaryEdges, and there is a point outside the loop which is either of type GAP_GRID
    // or it is an original surface point.  If the point to be used in interpolation is outside the loop,
    // then additional weighting is added equal to double the arc-length to be traversed going from the nearest
    // point on a fault to the end of a connected dying fault.
    // If setPoint is true, set the interpolated value and type at in,jn.
    // If mainDebug is true, compute the values but don't set them and print mainDebug output.
    public float interpolate(int in, int jn, boolean isInside, boolean setPoint, boolean debug)
    {
        byte pointType;
        byte fillFlag;
        boolean useGradient = true;
        boolean debugOutput;

        try
        {
            if(!insideGrid(in, jn)) return nullZValue;

            pointType = getPointType(in, jn);

            // debugOutput is true, if we want to mainDebug interpolation during construction
            debugOutput = doDebugInterpolate && in == debugInterpolateRow && jn == debugInterpolateCol;

            // debugInterpolate is true if we want mainDebug info, whether this is construction or subsequent debugging
            debugInterpolate = debug || debugOutput;

            // if mainDebug, we don't want to change existing values; we are redoing interpolation for mainDebug display
            if(!debug)
            {
                if(!isFillable(pointType)) return getPointZ(in, jn);

                if(isInside) fillFlag = GAP_FILLED;
                else fillFlag = GAP_EXTRAP;
            }
            else
            {
                debugInterpolatePointType = pointType;
                fillFlag = GAP_NONE;

                infoMessage("surface point mainDebug output. row: " + in +
                    " col: " + jn + " type: " + StsParameters.getGapTypeName(pointType));
            }

            float z = interpolator.interpolate(in, jn, useGradient, isInside, debugInterpolate);

            // setPoint is generally true, but if we are reading a surface position which is has not currently been
            // computed, we don't want to set values in z/type arrays

            z = setPointZ(in, jn, z, fillFlag);
            if(!debug && setPoint)
            {
                if(interpolator.getIsConverged())
                    setPointType(in, jn, fillFlag);
                else
                {
                    setPointType(in, jn, GAP_CANT_FILL);
                    interpolator.errorOutput();
                }
            }
            return z;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "leastSqFitInterpolate",
                " for row: " + in + " of " + pointRowMin + "-" + pointRowMax +
                " col: " + jn + " of " + pointColMin + "-" + pointColMax, e);

            setPointType(in, jn, GAP_CANT_FILL); // this sets blockGrid zValue to surface zValue
            return getPointZ(in, jn);
        }
    }

    public boolean delete()
    {
        super.delete();
        if(edgeLoop != null)
        {
            edgeLoop.delete();
            edgeLoop = null;
        }
        return true;
    }

    public void deleteTransientArrays()
    {
    }

    public void adjustThickness(StsBlockGrid baseBlockGrid, StsModelSurface baseSurface)
    {
        float topZ, baseZ;

        thickness = new float[nRows][nCols];

        for(int i = 0; i < nRows; i++)
        {
            for(int j = 0; j < nCols; j++)
            {
                topZ = getGapPointZ(i, j);
                if(topZ == nullZValue) continue;

                baseZ = baseBlockGrid.getGapPointZ(i, j);
                if(baseZ == nullZValue)
                {
                    baseZ = baseSurface.getZNotNull(i, j);
                    if(baseZ == nullZValue) continue;
                }
                if(topZ > baseZ)
                {
                    setGapPointZ(i, j, baseZ);
                    thickness[i][j] = 0.0f;
                }
                else
                    thickness[i][j] = baseZ - topZ;
            }
        }
    }

    public void interpolateThickness(StsModelSurface baseSurface, StsModelSurface topSurface,
                                     StsBlockGrid baseBlockGrid, StsBlockGrid topBlockGrid, float f)
    {
        float topZ, baseZ, z;

        for(int i = 0; i < nRows; i++)
        {
            for(int j = 0; j < nCols; j++)
            {
                topZ = topBlockGrid.getGapPointZ(i, j);
                if(topZ == nullZValue)
                {
                    topZ = topSurface.getZNotNull(i, j);
                    if(topZ == nullZValue) continue;
                }
                baseZ = baseBlockGrid.getGapPointZ(i, j);
                if(baseZ == nullZValue)
                {
                    baseZ = baseSurface.getZNotNull(i, j);
                    if(baseZ == nullZValue) continue;
                }

                if(topZ > baseZ)
                    setGapPointZ(i, j, baseZ);
                else
                    setGapPointZ(i, j, baseZ - f*(baseZ - topZ));
            }
        }
    }

    public byte getCellType(int row, int col)
    {
        byte cellType;
        byte[] cornerTypes = new byte[4];

        try
        {
            cornerTypes[0] = getPointType(row, col);
            cornerTypes[1] = getPointType(row+1, col);
            cornerTypes[2] = getPointType(row+1, col+1);
            cornerTypes[3] = getPointType(row, col+1);

            int nGapPnts = 0;
            int nGridPnts = 0;
            int nOutsidePnts = 0;

            for(int n = 0; n < 4; n++)
            {
                switch(cornerTypes[n])
                {
                    case GAP_GRID:                  // gridPoint for interpolating gapPoints
                        nGridPnts++;
                        break;
                    case GAP_NOT_FILLED:            // gapPoint we need to fill
                    case GAP_EXTRAP:                // extrapolated by least-sqs interpolated value
                    case GAP_FILLED:                // filled by least-sqs interpolated value
                    case GAP_FILL_HOLE:             // back-filled hole between edge and first GAP_GRID point
                    case GAP_CANT_FILL:             // couldn't fill (insufficient points?); assigned surface zValue
                        nGapPnts++;
                        break;
                    case GAP_NULL:                  // Not initialized
                    case GAP_SURF_GRID:             // outside domain
                        nOutsidePnts++;
                        break;
                    default:
                        StsException.systemError("StsBlockGrid.assignCellType() failed. " +
                        this.getLabel() + " row: " + row + " col: " + col +
                        " corner: " + n + " has undefined cellType: " + cornerTypes[n]);
                }
            }

            if(nGridPnts == 4)
                return CELL_FILL;
            else if(nOutsidePnts > 0)
                return CELL_EDGE;
            else
                return CELL_EMPTY;
        }
        catch(Exception e)
        {
            StsException.outputException("StsBlockGrid.assignCellType() failed.",
                e, StsException.WARNING);
            return CELL_EMPTY;
        }
    }

  	/** Not used currently: needs to be moved to StsLineLinkGrid */
/*
    public void drawSurfaceCursorEdge(StsGLPanel3d glPanel3d, StsXYSurfaceGridable grid,
                                             StsPoint point0, StsPoint point1, boolean drawGaps)
    {
	    float ftemp;
	    int itemp;
	    float dfx, dfy;
	    float slope;
	    int ix, iy;
	    float fx, fy;
	    int i, j;
        final int lineStart = 0;
        final int lineDraw = 1;
        final int lineNull = 2;
        int lineStatus = lineStart;

        StsGridPoint p0 = new StsGridPoint(point0, this);
        StsGridPoint p1 = new StsGridPoint(point1, this);
    	if(!currentModel.getGridDefinition().clipLineToGrid(p0, p1, this, false)) return;

        GL gl = glPanel3d.getGL();
        if(gl == null) return;

        gl.glLineWidth(StsGraphicParameters.edgeLineWidth);

	    int nColsm1 = nCols - 1;
	    int nRowsm1 = nRows - 1;

        StsModelSurfaceClass surfaceClass = (StsModelSurfaceClass)currentModel.getCreateStsClass(StsModelSurface.class);
        StsColor drawColor = surfaceClass.getDisplayHorizons() ? StsColor.BLACK : surface.getStsColor();
        drawColor.setGLColor(gl);
        StsColor nullColor = StsColor.WHITE;

    //	set_line_smoothing(current_window, TRUE);

	    StsPoint point = new StsPoint();
	    StsPoint lastPoint = new StsPoint();

	    if(p0.col == p1.col)		// map edge is along a col
	    {
		    if(p0.rowF > p1.rowF)
		    {
			    ftemp = p0.rowF;
			    p0.rowF = p1.rowF;
			    p1.rowF = ftemp;

			    itemp = p0.row;
			    p0.row = p1.row;
			    p1.row = itemp;
		    }

		    p0.row = (int)p0.rowF;
		    p1.row = (int)p1.rowF;

		    if(p0.col+1 <= nColsm1)
			    dfx = p0.colF - p0.col;
	    	else
		    {
			    p0.col = nColsm1 - 1;
			    dfx = 1.0f;
		    }

		    for(i = p0.row; i < p1.row; i++)
		    {
                if (rowInterpolateGridPoint(grid, i, p0.col, dfx, point))
                {
                    if (lineStatus == lineNull)  // draw the fault gap
                    {
                        lineStatus = lineStart;
                        if(drawGaps)
                        {
                            nullColor.setGLColor(gl);
                            gl.glBegin(GL.GL_LINE_STRIP);
                            gl.glVertex3fv(lastPoint.v, 0);
                            gl.glVertex3fv(point.v, 0);
                            gl.glEnd();
                            drawColor.setGLColor(gl);
                        }
                    }
                    if (lineStatus == lineStart) // start of regular line
                    {
                        gl.glBegin(GL.GL_LINE_STRIP);
                    }

                    gl.glVertex3fv(point.v, 0);  // draw regular line start or segment
                    lineStatus = lineDraw;
                    for (int p=0; p<3; p++) lastPoint.v[p] = point.v[p]; // save the point
                }
                else if (lineStatus == lineDraw)  // start of fault gap
                {
                    gl.glEnd();
                    lineStatus = lineNull;
                }
            }
            if (lineStatus == lineDraw) gl.glEnd(); // wrap up the last regular segment
        }
	    else if(p0.row == p1.row)		// map edge is along a row
	    {
		    if(p0.col > p1.col)
		    {
			    ftemp = p0.colF;
			    p0.colF = p1.colF;
			    p1.colF = ftemp;

			    itemp = p0.col;
			    p0.col = p1.col;
			    p1.col = itemp;
		    }

		    if(p0.row+1 <= nRowsm1)
			    dfy = p0.rowF - p0.row;
		    else
		    {
			    p0.row = nRowsm1 - 1;
			    dfy = 1.0f;
		    }

            gl.glBegin(GL.GL_LINE_STRIP);
		    for(j = p0.col; j < p1.col; j++)
		    {
                if (colInterpolateGridPoint(grid, p0.row, j, dfy, point))
                {
                    if (lineStatus == lineNull)  // draw the fault gap
                    {
                        lineStatus = lineStart;
                        if(drawGaps)
                        {
                            nullColor.setGLColor(gl);
                            gl.glBegin(GL.GL_LINE_STRIP);
                            gl.glVertex3fv(lastPoint.v, 0);
                            gl.glVertex3fv(point.v, 0);
                            gl.glEnd();
                            drawColor.setGLColor(gl);
                        }
                    }
                    if (lineStatus == lineStart) // start of regular line
                    {
                        gl.glBegin(GL.GL_LINE_STRIP);
                    }

                    gl.glVertex3fv(point.v, 0);  // draw regular line start or segment
                    lineStatus = lineDraw;
                    for (int p=0; p<3; p++) lastPoint.v[p] = point.v[p]; // save the point
                }
                else if (lineStatus == lineDraw)  // start of fault gap
                {
                    gl.glEnd();
                    lineStatus = lineNull;
                }
            }
            if (lineStatus == lineDraw) gl.glEnd(); // wrap up the last regular segment
	    }
	    else
	    {
		    int nc= Math.abs(p1.col - p0.col); // number of columns spanned
		    int nr = Math.abs(p1.row - p0.row); // number of rows spanned

		    if(nc > nr)
	    	{
			    if(p0.col > p1.col)
			    {
				    ftemp = p0.colF;
			    	p0.colF = p1.colF;
			    	p1.colF = ftemp;

				    itemp = p0.col;
			    	p0.col = p1.col;
			    	p1.col = itemp;

			    	ftemp = p0.rowF;
			    	p0.rowF = p1.rowF;
			    	p1.rowF = ftemp;

			    	itemp = p0.row;
			    	p0.row = p1.row;
			    	p1.row = itemp;
			    }

			    dfx = p0.colF - p0.col;
			    dfy = p0.rowF - p0.row;

			    slope = (p1.rowF - p0.rowF)/(p1.colF - p0.colF);

			    ix = p0.col + 1;

			    fy = p0.rowF + slope*(p0.col + 1 - p0.colF);
			    iy = (int)fy;

			    dfy = fy - iy;

                gl.glBegin(GL.GL_LINE_STRIP);
			    for(i = ix; i < p1.col; i++)
			    {
                    if (colInterpolateGridPoint(grid, iy, i, dfy, point))
                    {
                        if (lineStatus == lineNull)  // draw the fault gap
                        {
                            lineStatus = lineStart;
                            if(drawGaps)
                            {
                                nullColor.setGLColor(gl);
                                gl.glBegin(GL.GL_LINE_STRIP);
                                gl.glVertex3fv(lastPoint.v, 0);
                                gl.glVertex3fv(point.v, 0);
                                gl.glEnd();
                                drawColor.setGLColor(gl);
                            }
                        }
                        if (lineStatus == lineStart) // start of regular line
                        {
                            gl.glBegin(GL.GL_LINE_STRIP);
                        }

                        gl.glVertex3fv(point.v, 0);  // draw regular line start or segment
                        lineStatus = lineDraw;
                        for (int p=0; p<3; p++) lastPoint.v[p] = point.v[p]; // save the point
                    }
                    else if (lineStatus == lineDraw)  // start of fault gap
                    {
                        gl.glEnd();
                        lineStatus = lineNull;
                    }
                    fy += slope;
                    iy = (int)fy;
                    dfy = fy - iy;
                }
                if (lineStatus == lineDraw) gl.glEnd(); // wrap up the last regular segment
		    }
		    else
		    {
		    	if(p0.row > p1.row)
			    {
				    ftemp = p0.colF;
				    p0.colF = p1.colF;
				    p1.colF = ftemp;

				    itemp = p0.col;
				    p0.col = p1.col;
				    p1.col = itemp;

				    ftemp = p0.rowF;
				    p0.rowF = p1.rowF;
				    p1.rowF = ftemp;

				    itemp = p0.row;
				    p0.row = p1.row;
				    p1.row = itemp;
			    }

			    dfx = p0.colF - p0.col;
			    dfy = p0.rowF - p0.row;

			    slope = (p1.colF - p0.colF)/(p1.rowF - p0.rowF);

			    iy = p0.row + 1;

			    fx = p0.colF + slope*(p0.row + 1 - p0.rowF);
			    ix = (int)fx;

			    dfx = fx - ix;

                gl.glBegin(GL.GL_LINE_STRIP);
			    for(i = iy; i< p1.row; i++)
			    {
                    if (rowInterpolateGridPoint(grid, i, ix, dfx, point))
                    {
                        if (lineStatus == lineNull)  // draw the fault gap
                        {
                            lineStatus = lineStart;
                            if(drawGaps)
                            {
                                nullColor.setGLColor(gl);
                                gl.glBegin(GL.GL_LINE_STRIP);
                                gl.glVertex3fv(lastPoint.v, 0);
                                gl.glVertex3fv(point.v, 0);
                                gl.glEnd();
                                drawColor.setGLColor(gl);
                            }
                        }
                        if (lineStatus == lineStart) // start of regular line
                        {
                            gl.glBegin(GL.GL_LINE_STRIP);
                        }

                        gl.glVertex3fv(point.v, 0);  // draw regular line start or segment
                        lineStatus = lineDraw;
                        for (int p=0; p<3; p++) lastPoint.v[p] = point.v[p]; // save the point
                    }
                    else if (lineStatus == lineDraw)  // start of fault gap
                    {
                        gl.glEnd();
                        lineStatus = lineNull;
                    }
				    fx += slope;
				    ix = (int)fx;
				    dfx = fx - ix;
                }
                if (lineStatus == lineDraw) gl.glEnd(); // wrap up the last regular segment
		    }
	    }
    //	set_line_smoothing(current_window, FALSE);
    }

    private boolean rowInterpolateGridPoint(StsXYSurfaceGridable grid, int i, int j, float f, StsPoint point)
    {
        if (point==null) return false;

        StsGridPoint point0 = new StsGridPoint(i, j, grid);
        StsGridPoint point1 = new StsGridPoint(i+1, j, grid);

        if(f <= 0.0f)
        {
            if(point0 == null || !point0.notNull()) return false;
            point.setValues(point0);
            return true;
        }
        else if(f >= 1.0f)
        {
            if(point1 == null || !point1.notNull()) return false;
            point.setValues(point1);
            return true;
        }
        else
        {
            if(point0 == null || !point0.notNull()) return false;
            if(point1 == null || !point1.notNull()) return false;

            point.interpolatePoints(point0.pxyz, point1.pxyz, f);
            return true;
        }
    }

    private boolean colInterpolateGridPoint(StsXYSurfaceGridable grid, int i, int j, float f, StsPoint point)
    {
        if (point==null) return false;

        StsGridPoint point0 = new StsGridPoint(i, j, grid);
        StsGridPoint point1 = new StsGridPoint(i, j+1, grid);

        if(f <= 0.0f)
        {
            if(point0 == null || !point0.notNull()) return false;
            point.setValues(point0);
            return true;
        }
        else if(f >= 1.0f)
        {
            if(point1 == null || !point1.notNull()) return false;
            point.setValues(point1);
            return true;
        }
        else
        {
            if(point0 == null || !point0.notNull()) return false;
            if(point1 == null || !point1.notNull()) return false;

            point.interpolatePoints(point0, point1, f);
            return true;
        }
   }
*/
    public void getUnboundedGridFloatCoordinates(StsPoint point, StsPoint2D ijPoint)
    {
        getUnboundedGridFloatCoordinates(point.v, ijPoint);
    }

    private void getUnboundedGridFloatCoordinates(StsPoint2D point, StsPoint2D ijPoint)
    {
        getUnboundedGridFloatCoordinates(point, ijPoint);
    }

    public StsPoint2D getUnboundedGridFloatCoordinates(StsPoint point)
    {
        StsPoint2D ijPoint = new StsPoint2D();
        getUnboundedGridFloatCoordinates(point.v, ijPoint);
        return ijPoint;
    }

    private void getUnboundedGridFloatCoordinates(float[] xy, StsPoint2D ijPoint)
    {
        ijPoint.x = getColCoor(xy);
        ijPoint.y = getRowCoor(xy);

    }

    /** find the surface point at a mouse cursor location */
    public StsGridSectionPoint getSurfacePosition(StsPoint gridPoint0, StsPoint gridPoint1)
    {
        boolean debug = false; //StsTrace.getTrace();

        StsPoint point0 = new StsPoint(gridPoint0);
        StsPoint point1 = new StsPoint(gridPoint1);
/*
        float xMin = surface.getColX(pointColMin);
        float xMax = surface.getColX(pointColMax);
        float yMin = surface.getRowY(pointRowMin);
        float yMax = surface.getRowY(pointRowMax);

        if(!StsCube.clipLineToCube(point0.v, point1.v, xMin, xMax, yMin, yMax, zMin, zMax)) return null;
*/
        // calculate a increment to step along the near-far line
        float df = 0.5f / Math.max(nCols, nRows);

        // step along the near-far line looking for an intersection
        float lastError = nullZValue;
        StsPoint point = new StsPoint(5);
        StsPoint lastPoint = new StsPoint(5);

        float f = 0.0f;
        while (f < 1.0f)
        {
            point.interpolatePoints(point0, point1, f);
            if (debug) System.out.println("interp Pt:  x = " + point.v[0] +
                ", y = " + point.v[1] + ", z = " + point.v[2]);
            float z = interpolateBilinearZ(point, false, false);
            if (debug) System.out.println("interpolated z = " + z);
            if (z != nullZValue)
            {
                float error = z - point.v[2];
                if (error==0.0f) break;
                if (lastError != nullZValue)
                {
                    if (lastError<0.0f && error>0.0f ||
                        lastError>0.0f && error<0.0f)
                    {
                        f = Math.abs(lastError/(lastError-error));
                        point0 = lastPoint;
                        point1 = point;
                        point.interpolatePoints(point0, point1, f);
                        if (debug) System.out.println("last interp Pt:  x = " + point.v[0] +
                            ", y = " + point.v[1] + ", z = " + point.v[2]);

                        return getSurfacePosition(point);
                    }
                }
                lastError = error;
                for (int i=0; i<3; i++) lastPoint.v[i] = point.v[i];
            }
            f += df;
        }

        return null;
    }

    public StsGridSectionPoint getSurfacePosition(StsPoint p)
    {
        // Get ij location info
        StsPoint2D ij = new StsPoint2D();
        getUnboundedGridFloatCoordinates(p.v, ij);

        // display x-y-z location

        int i = (int)ij.getY();
        int j = (int)ij.getX();
        byte pointType = getPointType(i, j);

        if(isPointNotActive(pointType)) return null;
        StsGridSectionPoint gridPoint = new StsGridSectionPoint(p, i, j, null, this, false);
//        gridPoint.setEdgeGroup(edgeGroup);
        return gridPoint;
    }

    public void displaySurfacePosition(StsGridPoint gridPoint)
    {
        int i = gridPoint.row;
        int j = gridPoint.col;
        byte pointType = getPointType(i, j);
        float[] xyz = gridPoint.getXYZorT();

        String message = null;
        if (surface.getPropertyVolume() == null)
        {
            message = new String("Surface: " + surface.getName() +
                                 "Block grid: " + getIndex() +
                                 "   row: " + i + " col: " + j +
                                 "  X: " + (int)xyz[0] +
                                 "  Y: " + (int)xyz[1] +
                                 "  Z: " + (int)xyz[2] +
                                 " type: " + StsParameters.getGapTypeName(pointType) );
        }
        else
        {
            String propertyName = surface.getPropertyName();
            String property = propertyName==null ? "  Property: " : "  " + propertyName+": ";
            float value = surface.getPropValue(xyz);
            message = new String("Surface: " + surface.getName() +
                                 "Block grid: " + getIndex() +
                                 "   row: " + i + " col: " + j +
                                 "  X: " + (int)xyz[0] +
                                 "  Y: " + (int)xyz[1] +
                                 "  Z: " + (int)xyz[2] +
                                 property + value +
                                 " type: " + StsParameters.getGapTypeName(pointType) );
        }
        StsMessageFiles.infoMessage(message);
    }

    public void debugSurfacePoint(StsGridPoint gridPoint)
    {
        int row = gridPoint.row;
        int col = gridPoint.col;

        byte pointType = getPointType(row, col);

        switch(pointType)
        {
            case GAP_FILLED:
            case GAP_FILL_HOLE:
                debugSurfacePoint(row, col, true);
                break;
            case GAP_CUT:
            case GAP_EXTRAP:
            case GAP_CANT_FILL:
                debugSurfacePoint(row, col, false);
                break;
            default:
                debugSurfacePointMessage(row, col, pointType);
        }
    }

    private void debugSurfacePointMessage(int row, int col, byte pointType)
    {
        StsMessageFiles.infoMessage("surface point mainDebug output\n" + " row: " + row +
                    " col: " + col + " type: " + StsParameters.getGapTypeName(pointType) + "\n" +
                    " no additional info.");
    }

    private void debugSurfacePoint(int row, int col, boolean isInside)
    {
        debugInterpolateRow = row;
        debugInterpolateCol = col;
        // StsBlockGridInterpolationWeightedPlane.getInstance(surface, this);
        interpolate(row, col, isInside, SET_POINT_FALSE, DEBUG_TRUE);
    }

    private String getTypeName(int row, int col)
    {
        return StsParameters.getGapTypeName(getPointType(row, col));
    }

    public float getPointZ(int row, int col)
    {
        if(pointsZ == null) return nullZValue;
        if(insideGrid(row, col))
        {
            float z = pointsZ[row-pointRowMin][col-pointColMin];
            if(z != nullValue) return z;
        }
        if(surface.insideGrid(row, col)) return surface.getZorT(row, col);

        StsException.systemError("StsBlockGrid.getPointZ() failed." + outsideGridMessage(row, col));
        return nullZValue;
    }

    public StsPoint getComputePoint(float rowF, float colF)
    {
        StsGridPoint gridPoint = new StsGridPoint(rowF, colF, this);
        if(!interpolateBilinear(gridPoint, true, true)) return null;
        return new StsPoint(gridPoint.getPoint());
    }

    public float getComputePointZ(float rowF, float colF)
    {
        StsGridPoint gridPoint = new StsGridPoint(rowF, colF, this);
        if (!interpolateBilinear(gridPoint, true, true)) return StsParameters.nullValue;
        return gridPoint.getPoint().getZorT();
    }

    public float getPointZ(float rowF, float colF)
    {
        StsGridPoint gridPoint = new StsGridPoint(rowF, colF, this);
        if (!interpolateBilinear(gridPoint, false, false)) return StsParameters.nullValue;
        return gridPoint.getPoint().getZorT();
    }

    public StsPoint getComputePoint(int row, int col, boolean isInside)
    {
        float x, y, z;
        z = interpolate(row, col, isInside, SET_POINT_TRUE, DEBUG_FALSE);
        if(z == nullZValue) return null;
        x = getXCoor(row, col);
        y = getYCoor(row, col);
        return new StsPoint( x, y, z );
    }

    public float getComputePointZ(int row, int col)
    {
        return interpolateOutsidePoint(row, col, true);
    }

    public float getComputeExtrapolatedZ(int row, int col)
    {
        StsBlockGridInterpolationWeightedPlane.extrapolator(surface, this);
        return interpolateOutsidePoint(row, col, true);
    }

    public float getComputePointZ(int row, int col, boolean computeIfNull, boolean setPoint)
    {
        if(!computeIfNull) return getPointZ(row, col);
        return interpolateOutsidePoint(row, col, setPoint);
     /*
        float z = getPointZ(row, col);
        if(z != nullZValue || !computeIfNull) return z;
        return interpolateOutsidePoint(row, col, setPoint);
     */
    }
/*
    public float getComputePointZ(int row, int col, boolean computeIfNull, boolean setPoint)
    {
        if(!computeIfNull) return getPointZ(row, col);
        float z = leastSqFitInterpolateOutsidePoint(row, col, setPoint); // z is time or approx-depth
        StsSeismicVelocityModel velocityModel = currentModel.project.getSeismicVelocityModel();
        if(velocityModel == null) return nullValue;
        float x = getXCoor((float)row, (float)col);
        float y = getYCoor((float)row, (float)col);
        float d = (float)velocityModel.getZ(x, y, z);
        this.adjPointsZ[row][col] = d;
        return z;
    }
*/
   	public float interpolateBilinearZ(StsPoint point, boolean computeIfNull, boolean setPoint)
	{
        if (point == null) return nullZValue;

		StsGridPoint gridPoint = new StsGridPoint(point, this);
        boolean OK = interpolateBilinear(gridPoint, computeIfNull, setPoint);
		if(OK)
            return gridPoint.getZ();
		else
			return nullValue;
	}

    public float interpolateBilinearZ(StsGridPoint gridPoint, boolean computeIfNull, boolean setPoint)
    {
        if(!interpolateBilinear(gridPoint, computeIfNull, setPoint)) return nullValue;
        return gridPoint.getZ();
    }

    private boolean interpolateBilinear(StsGridPoint gridPoint, boolean compute, boolean setPoint)
    {
		float rowF, colF;
        int i, j;
        float dx, dy;
        float z, w, weight = 0, zWeighted = 0;
        boolean allNulls = true;
        final boolean debug = false; //StsTrace.getTrace();

		rowF = gridPoint.rowF;
		colF = gridPoint.colF;
        if(!insideGrid(rowF, colF)) return false;

        i = (int)rowF;
        if(i == pointRowMax) i--;
        dy = rowF - i;

        j = (int)colF;
        if(j == pointColMax) j--;
        dx = colF - j;

        if (debug) System.out.println("\trow = " + rowF +
                    ", dy = " + dy + ", col = " + colF +
                    ", dx = " + dx);

        if((w = (1.0f-dy) * (1.0f-dx)) > StsParameters.roundOff &&
           (z = getComputePointZ(i, j, compute, setPoint)) != nullZValue)
        {
            weight += w;
            zWeighted += w*z;
            allNulls = false;
            if (debug) System.out.println("\tz[i][j] = " + z);
        }
        if((w = dy * (1.0f-dx)) > StsParameters.roundOff &&
           (z = getComputePointZ(i+1, j, compute, setPoint)) != nullZValue)
        {
            weight += w;
            zWeighted += w*z;
            allNulls = false;
            if (debug) System.out.println("\tz[i+1][j] = " + z);
        }
        if((w = (1.0f-dy) * dx) > StsParameters.roundOff &&
           (z = getComputePointZ(i, j+1, compute, setPoint)) != nullZValue)
        {
            weight += w;
            zWeighted += w*z;
            allNulls = false;
            if (debug) System.out.println("\tz[ijE] = " + z);
        }
        if((w = dy*dx) > StsParameters.roundOff &&
           (z = getComputePointZ(i+1, j+1, compute, setPoint)) != nullZValue)
        {
            weight += w;
            zWeighted += w*z;
            allNulls = false;
            if (debug) System.out.println("\tz[i+1][j+1] = " + z);
        }

        if (allNulls) return false;

        if(weight < 1.0f) zWeighted /= weight;
        gridPoint.setZ(zWeighted);
		return true;
    }

    public StsPoint getStsPoint(int row, int col)
    {
        float x, y, z;
        z = getPointZ(row, col);
        if(z == nullZValue) z = surface.getZorT(row, col);
        x = getXCoor(row, col);
        y = getYCoor(row, col);
        return new StsPoint(x, y, z);
    }

    private float setPointZ(int row, int col, float z, byte fillFlag)
    {
        boolean debugThisPoint = Main.debugPoint && row == Main.debugPointRow && col == Main.debugPointCol;

        if(debugThisPoint)
        {
            StsException.systemDebug("StsBlockGrid.setPointZ() called. " +
                this.getLabel() + " row: " + row + " col: " + col + " pointsType: " +
                StsParameters.getGapTypeName(pointsType[row-pointRowMin][col-pointColMin]) +
                " z set to: " + z);
        }

        if(pointsZ == null) return z;

        if(fillFlag != GAP_FILLED && fillFlag != GAP_EXTRAP)
            z = checkIntersectAgainstZBelow(row, col, z, debugThisPoint);
        pointsZ[row-pointRowMin][col-pointColMin] = z;
        return z;
    }

    public float checkIntersectAgainstZBelow(int row, int col, float z, boolean debugThisPoint)
    {
        if(blockGridBelow == null || blockGridBelow.getPointsZ() == null) return z;
        float zBelow = blockGridBelow.getPointZ(row, col);
        if(zBelow == nullValue || z < zBelow) return z;

        float adjustedZ = zBelow - 0.5f;
        if(debugThisPoint) StsException.systemDebug("StsBlockGrid.setPointZ() called. " +
                this.getLabel() + " row: " + row + " col: " + col + " pointsType: " +
                StsParameters.getGapTypeName(pointsType[row-pointRowMin][col-pointColMin]) +
                " adjusted from: " + z + " to: " + adjustedZ);
        return adjustedZ;
    }

    public final float[] getNormal(int row, int col)
    {
        if(insideGrid(row, col))
            return normals[row-pointRowMin][col-pointColMin];
        else
            return vertNormal;
    }

    public void checkConstructGridNormals()
    {
        if(normals != null) return;
        computeGridNormals();
    }

   	public void computeGridNormals()
    {
		int row = -1, col = -1;
        float z, zRite, zUp;
        float[] normal = null;

        int nNewRows = pointRowMax - pointRowMin + 1;
        int nNewCols = pointColMax - pointColMin + 1;
        normals = new float[nNewRows][nNewCols][3];

        float[] difRite = new float[3];
        float[]  difUp = new float[3];

        difUp[0] = 0.0f;
        difUp[1] = surface.getYInc();
        difRite[0] = surface.getXInc();
        difRite[1] = 0.0f;

	    try
	    {
			for(row = pointRowMin; row < pointRowMax; row++)
			{
				zRite = getPointZ(row, pointColMin);
				for(col = pointColMin; col < pointColMax; col++)
				{
					z = zRite;
					zRite = getPointZ(row, col + 1);
					if(zRite == nullZValue)zRite = surface.getZorT(row, col + 1);
					zUp = getPointZ(row + 1, col);
					if(zUp == nullZValue)zUp = surface.getZorT(row + 1, col);
					difUp[2] = zUp - z;
					difRite[2] = zRite - z;
					normal = StsMath.crossProduct(difUp, difRite);
					addNormal(row, col, normal);
//					setNormal(row, col, normal);
				}
				// Copy normal for last column from second to last column
				setNormal(row, pointColMax, normal);
			}

			// Copy normals for top row from second to top row
			for(col = pointColMin; col <= pointColMax; col++)
			{
				normal = getNormal(pointRowMax - 1, col);
				setNormal(pointRowMax, col, normal);
			}
		}
		catch(Exception e)
		{
			StsException.outputException("StsBlockGrid.computeGridNormals() failed at row " + row + " col " + col, e, StsException.WARNING);
		}
	}

	private void addNormal(int row, int col, float[] normal)
	{
		try
		{
			row -= pointRowMin;
			col -= pointColMin;
			StsMath.vectorAdd(normals[row][col], normal);
			StsMath.vectorAdd(normals[row + 1][col], normal);
			StsMath.vectorAdd(normals[row][col + 1], normal);
		}
		catch(Exception e)
		{
			StsException.outputException("StsBlockGrid.addNormal() failed for row " + row + " col " + col, e, StsException.WARNING);
		}
	}

    public void setNormal(int row, int col, float[] normal)
    {
        if(normal == null || normals == null) return;
		try
		{
			normals[row - pointRowMin][col - pointColMin] = normal;
		}
		catch(Exception e)
		{
			StsException.outputException("StsBlockGrid.setNormal() failed for row " + row + " col " + col, e, StsException.WARNING);
		}
	}

    public float[] getNormal(StsGridSectionPoint point)
    {
        return getNormal(point.getRowF(null), point.getColF(null));
    }

    public float[] getNormal(float rowF, float colF)
    {
        float[] norm0, norm1, norm2;

        try
        {
            int row = (int)rowF;
            float dR = rowF - row;

            int col = (int)colF;
            float dC = colF - col;

            if(dR == 0.0f)
            {
                if(dC == 0.0f)
                    return getNormal(row, col);
                else
                    return StsMath.interpolate(getNormal(row, col), getNormal(row, col+1), dC);
            }
            else
            {
                norm0 = StsMath.interpolate(getNormal(row, col), getNormal(row+1, col), dR);
                if(dC == 0.0f)
                    return norm0;
                else
                {
                    norm1 = StsMath.interpolate(getNormal(row, col+1), getNormal(row+1, col+1), dR);
                    return StsMath.interpolate(norm0, norm1, dC);
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsBlockGrid.getNormal() failed.", e, StsException.WARNING);
            return null;
        }
    }

    public float getGapPointZ(int row, int col)
    {
        if(isPointGap(row, col))
            return getPointZ(row, col);
        else
            return nullZValue;
    }

    public void setGapPointZ(int i, int j, float z)
    {
        setPointZ(i, j, z, GAP_NOT_FILLED);
        setPointType(i, j, GAP_NOT_FILLED);
    }

    /* get fractional map indices for a given x, y (or false) */
    private boolean getMapIndices(StsPoint p, StsPoint2D floatXYIndex)
    {
        getUnboundedGridFloatCoordinates(p, floatXYIndex);

        float rowF = floatXYIndex.y;
        float colF = floatXYIndex.x;

//        if(!insideGrid(rowF, colF)) return false;

        return true;
    }

    public void getNearestMapIndices(StsPoint p, StsPoint2D floatXYIndex)
    {
        getUnboundedGridFloatCoordinates(p, floatXYIndex);

        // adjust to nearest point
        floatXYIndex.x = StsMath.minMax(floatXYIndex.x, 0.0f, (float)nCols-1);
        floatXYIndex.y = StsMath.minMax(floatXYIndex.y, 0.0f, (float)nRows-1);
    }

    // downsize the blockGrids from project grid size to bounding box with padding of 2 on each side;
    // this might not be still working with the subsized starting array
    public void downSizePointArrays()
    {
        int row, col;

        int newRowMin = Math.max(newPointRowMin-1, 0);
        int newRowMax = Math.min(newPointRowMax+1, nRows-1);
        int newColMin = Math.max(newPointColMin-1, 0);
        int newColMax = Math.min(newPointColMax+1, nCols-1);

        int nNewRows = newRowMax-newRowMin+1;
        int nNewCols = newColMax-newColMin+1;

        float[][] newPointsZ = new float[nNewRows][nNewCols];
        for(row = newRowMin; row <= newRowMax; row++)
        {
            float[] arrayRow = pointsZ[row-pointRowMin];
            float[] newArrayRow = newPointsZ[row-newRowMin];
            System.arraycopy(arrayRow, newColMin, newArrayRow, 0, nNewCols);
        }
        pointsZ = newPointsZ;

        byte[][] newPointsType = new byte[nNewRows][nNewCols];
        for(row = newRowMin; row <= newRowMax; row++)
        {
            byte[] arrayRow = pointsType[row-pointRowMin];
            byte[] newArrayRow = newPointsType[row-newRowMin];
            System.arraycopy(arrayRow, newColMin, newArrayRow, 0, nNewCols);
        }
        pointsType = newPointsType;

        pointRowMin = newRowMin;
        pointRowMax = newRowMax;
        pointColMin = newColMin;
        pointColMax = newColMax;

        nRows = pointRowMax - pointRowMin + 1;
        nCols = pointColMax - pointColMin + 1;
    }

    // blockGrid contains points which extend beyond domain: these are bounded by row/col min/max;
    // domainGrid limits are bounded by domain row/col min/max;
/*
    private void upSizePointArrays()
    {
        int nNewRows, nNewCols;
        int row, col;

        // resize blockGrid arrays to bounds of blockGrid points and domain points

        nNewRows = newPointRowMax-newPointRowMin+1;
        nNewCols = newPointColMax-newPointColMin+1;

        int nRowPnts = newPointColMax-newPointColMin+1;
        int oldColStart = pointColMin-newPointColMin;

        float[][] newPointsZ = new float[nNewRows][nNewCols];

        for(row = 0; row < nNewRows; row++)
            for(col = 0; col < nNewCols; col++)
                newPointsZ[row][col] = nullZValue;

        for(row = newPointRowMin; row <= newPointRowMax; row++)
        {
            float[] arrayRow = pointsZ[row-pointRowMin];
            float[] newArrayRow = newPointsZ[row-newPointRowMin];
            System.arraycopy(arrayRow, oldColStart, newArrayRow, 0, nRowPnts);
        }
        pointsZ = newPointsZ;

        byte[][] newPointsType = new byte[nNewRows][nNewCols];
        for(row = newPointRowMin; row <= newPointRowMax; row++)
        {
            byte[] arrayRow = pointsType[row-pointRowMin];
            byte[] newArrayRow = newPointsType[row-newPointRowMin];
            System.arraycopy(arrayRow, oldColStart, newArrayRow, 0, nRowPnts);
        }
        pointsType = newPointsType;

        checkSetSurfaceNonNulls();

        // resize domain cell limits

        pointRowMin = newPointRowMin;
        pointRowMax = newPointRowMax;
        pointColMin = newPointColMin;
        pointColMax = newPointColMax;

        nRows = nNewRows;
        nCols = nNewCols;
    }
*/
/*
    private void checkSetSurfaceNonNulls()
    {
        for(int row = pointRowMin; row <= pointRowMax; row++)
            for(int col = pointColMin; col <= pointColMax; col++)
                if(isPointType(row, col, GAP_NULL) && !surface.isPointNull(row, col))
                    setPointType(row, col, GAP_SURF_GRID);
    }
*/
    // if we are on a row we want the pointColMin and versa visa
    public int getCrossingRowColMin(int rowOrCol)
    {
        if(rowOrCol == ROW)
            return pointColMin;
        else if(rowOrCol == COL)
            return pointRowMin;
        else
            return -1;
    }

    public int getCrossingRowColMax(int rowOrCol)
    {
        if(rowOrCol == ROW)
            return pointColMax;
        else if(rowOrCol == COL)
            return pointRowMax;
        else
            return -1;
    }

    /** find the surface point at a mouse cursor location */
    public StsPoint getPointOnSurface(StsPoint gridPoint0, StsPoint gridPoint1)
    {
        boolean debug = false; //StsTrace.getTrace();

        StsPoint point0 = new StsPoint(gridPoint0);
        StsPoint point1 = new StsPoint(gridPoint1);

//        if(!StsCube.clipLineToCube(point0.v, point1.v, xMin, xMax, yMin, yMax, zMin, zMax)) return null;

        // calculate a increment to step along the near-far line
        float df = 0.5f / Math.max(nCols, nRows);

        // step along the near-far line looking for an intersection
        float lastError = nullZValue;
        StsPoint point = new StsPoint(5);
        StsPoint lastPoint = new StsPoint(5);

        float f = 0.0f;
        while (f < 1.0f)
        {
            point.interpolatePoints(point0, point1, f);
            if (debug) System.out.println("interp Pt:  x = " + point.v[0] +
                ", y = " + point.v[1] + ", z = " + point.v[2]);
            float z = interpolateBilinearZ(point, false, false);
            if (debug) System.out.println("interpolated z = " + z);
            if (z != nullZValue)
            {
                float error = z - point.v[2];
                if (error==0.0f) break;
                if (lastError != nullZValue)
                {
                    if (lastError<0.0f && error>0.0f ||
                        lastError>0.0f && error<0.0f)
                    {
                        f = Math.abs(lastError/(lastError-error));
                        point0 = lastPoint;
                        point1 = point;
                        point.interpolatePoints(point0, point1, f);
                        if (debug) System.out.println("last interp Pt:  x = " + point.v[0] +
                            ", y = " + point.v[1] + ", z = " + point.v[2]);

                        // final range checks - do we need to restore these?? TJL 10/21/02
                    /*
                        if (point.v[0]<xMin) point.v[0] = xMin;
                        else if (point.v[0]>xMax) point.v[0] = xMax;
                        if (point.v[1]<yMin) point.v[1] = yMin;
                        else if (point.v[1]>yMax) point.v[1] = yMax;
                    */
 //                       if (point.v[2]<zMin) point.v[2] = zGapMin;
 //                       else if (point.v[2]>zMax) point.v[2] = zMax;
                        return point;
                    }

                }
                lastError = error;
                for (int i=0; i<3; i++) lastPoint.v[i] = point.v[i];
            }
            f += df;
        }

        return null;
    }

    public void adjustAuxiliaryEdges()
    {
        StsList edges = getEdges();
        int nEdges = edges.getSize();
        for(int e = 0; e < nEdges; e++)
        {
            StsSurfaceEdge edge = (StsSurfaceEdge)edges.getElement(e);
            if(edge.isAuxiliary()) adjustAuxiliaryEdges(edge);
        }
    }

    private void adjustAuxiliaryEdges(StsSurfaceEdge edge)
    {
        StsEdgeLoopRadialGridLink minusLink, plusLink;

        Iterator<StsEdgeLoopRadialGridLink> edgeLinksIterator = edgeLoop.getEdgeLinksIterator(edge, PLUS);
        while(edgeLinksIterator.hasNext())
        {
            minusLink = edgeLinksIterator.next();
            int outsideLinkIndex = getOutsidePlusLinkIndex(minusLink);
            if(outsideLinkIndex == StsParameters.NO_LINK) continue;
            plusLink = minusLink.getLink(outsideLinkIndex);
            int rowOrCol = StsEdgeLoopRadialGridLink.getLinkIndexRowOrCol(outsideLinkIndex);
            adjustAuxiliaryEdges(rowOrCol, minusLink, plusLink);
        }
    }

    private void adjustAuxiliaryEdges(int rowOrCol, StsEdgeLoopRadialGridLink minusLink,
                                                         StsEdgeLoopRadialGridLink plusLink)
    {
        StsGridSectionPoint minusPoint = null, plusPoint = null;
        StsBlockGrid minusBlockGrid, plusBlockGrid;

        if(plusLink == null) return;

        try
        {
			minusPoint = minusLink.getPoint();
			plusPoint = plusLink.getPoint();

            minusBlockGrid = (StsBlockGrid)minusLink.getEdgeLoop().getGrid();
            plusBlockGrid = (StsBlockGrid)plusLink.getEdgeLoop().getGrid();

            int minusRow = (int)minusLink.getRowColF(ROW);
            int minusCol = (int)minusLink.getRowColF(COL);
            int plusRow, plusCol;

            if(rowOrCol == ROW)
            {
                plusRow = minusRow;
                plusCol = minusCol+1;
            }
            else
            {
                plusRow = minusRow+1;
                plusCol = minusCol;
            }

            if(minusBlockGrid.isPointType(minusRow, minusCol, GAP_NULL)) return;
            if(plusBlockGrid.isPointType(plusRow, plusCol, GAP_NULL)) return;

            float zMinus = minusBlockGrid.getPointZ(minusRow, minusCol);
            float zPlus = plusBlockGrid.getPointZ(plusRow, plusCol);

            float rowColF = minusPoint.getCrossingRowColF(null, rowOrCol);
            float f = rowColF - (int)rowColF;

            float z = zMinus + f*(zPlus - zMinus);
            minusPoint.getPoint().setZ(z);
            plusPoint.getPoint().setZ(z);
        }
        catch(Exception e)
        {
            StsException.outputException("StsEdgeGroup.adjustAuxiliaryEdges() failed." +
                " for minusPoint: " + StsGridSectionPoint.staticGetLabel(minusPoint) +
                " and plusPoint: " +  StsGridSectionPoint.staticGetLabel(plusPoint),
                e, StsException.WARNING);
        }
    }

    public void adjustAuxiliaryEdges(StsBlockGrid minusBlockGrid, int minusRow, int minusCol,
                    StsBlockGrid plusBlockGrid, int plusRow, int plusCol)
    {
        byte pointType;
        float z;

        try
        {
            if(plusBlockGrid.insideGrid(minusRow, minusCol))
            {
                pointType = minusBlockGrid.getPointType(minusRow, minusCol);
                if(pointType != GAP_NULL)
                {
                    z = minusBlockGrid.getPointZ(minusRow, minusCol);
                    plusBlockGrid.setPointType(minusRow, minusCol, pointType);
                    plusBlockGrid.setPointZ(minusRow, minusCol, z, pointType);
                }
            }
            if(minusBlockGrid.insideGrid(plusRow, plusCol))
            {
                pointType = plusBlockGrid.getPointType(plusRow, plusCol);
                if(pointType != GAP_NULL)
                {
                    z = plusBlockGrid.getPointZ(plusRow, plusCol);
                    minusBlockGrid.setPointType(plusRow, plusCol, pointType);
                    minusBlockGrid.setPointZ(plusRow, plusCol, z, pointType);
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsEdgeGroup.adjustAuxiliaryEdges() failed.",
                e, StsException.WARNING);
        }
    }

    private int getOutsidePlusLinkIndex(StsEdgeLoopRadialGridLink minusLink)
    {
        StsEdgeLoopRadialGridLink plusLink;

        int rowOrCol = minusLink.getRowOrCol();

        if( (rowOrCol == ROW || rowOrCol == ROWCOL) &&
            minusLink.hasInsideDirection(ROW, MINUS) && minusLink.hasOutsideDirection(ROW, PLUS))
                return StsEdgeLoopRadialGridLink.ROW_PLUS;

        if((rowOrCol == COL || rowOrCol == ROWCOL) &&
            minusLink.hasInsideDirection(COL, MINUS) && minusLink.hasOutsideDirection(COL, PLUS))
                return StsEdgeLoopRadialGridLink.COL_PLUS;

        return StsEdgeLoopRadialGridLink.NO_LINK;
    }

    // this routine needs to be fixed now that cellTypes aren't available
    public boolean fillZoneBlockArray(StsZoneBlock[][] zoneBlockArray,
                        int rowMin, int rowMax, int colMin, int colMax,
                        StsZoneBlock zoneBlock)
    {
        try
        {
            for(int i = rowMin; i < rowMax; i++)
                for(int j = colMin; j < colMax; j++)
                      zoneBlockArray[i-rowMin][j-colMin] = zoneBlock;
//                    if(isCellTypeQuad(i, j)) zoneBlockArray[i-rowMin][j-colMin] = zoneBlock;
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsBlockGrid.fillZoneBlockArray() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    public void fillInsideGridPoints()
    {
		StsEdgeLoopRadialGridLink link;
        // do this with linkGrid instead...
        try
        {
            StsList loopLinks = edgeLoop.getLoopLinks();
            int nLinks = loopLinks.getSize();

            for(int n = 0; n < nLinks; n++)
            {
                link = (StsEdgeLoopRadialGridLink)loopLinks.getElement(n);
                if(link.isConnected(ROW, PLUS)) fillInsideGridPoints(link, ROW);
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsBlockGrid.fillInsidePoints() failed.",
                e, StsException.WARNING);
        }
    }

    private void fillInsideGridPoints(StsEdgeLoopRadialGridLink link, int rowOrCol)
    {
        int rowColStart, rowColEnd;
        int nPnts;

        if(link == null) return;
        StsEdgeLoopRadialGridLink nextLink = link.getLink(rowOrCol, PLUS);
        if(nextLink == null) return;
        if(!nextLink.isConnected(rowOrCol, MINUS)) return;  // should be true

        rowColStart = StsMath.ceiling(link.getCrossingRowColF(rowOrCol));
        rowColEnd = StsMath.floor(nextLink.getCrossingRowColF(rowOrCol));
        nPnts = rowColEnd - rowColStart + 1;

        if(nPnts <= 0) return;

        int rowInc, colInc;
        int row, col;

        if(rowOrCol == ROW)
        {
            rowInc = 0;
            colInc = 1;
            row = link.getRowCol(ROW);
            col = rowColStart;
        }
        else
        {
            rowInc = 1;
            colInc = 0;
            row = rowColStart;
            col = link.getRowCol(COL);
        }

        for(int n = 0; n < nPnts; n++)
        {
			if(!insideGrid(row, col)) continue;

			byte pointType = getPointType(row, col);
            if(pointType == GAP_NULL && surface.isPointNotNull(row, col))
                setPointType(row, col, GAP_GRID);
             else if(pointType == GAP_NOT_FILLED)
                interpolate(row, col, true, true, false);

        /*
            if(pointType == GAP_NULL)
            {
                if(surface.isPointNull(row, col))
                    leastSqFitInterpolate(row, col, true, true, false);
                else
                    setPointType(row, col, GAP_GRID);
             }
             else if(pointType == GAP_NOT_FILLED)
                leastSqFitInterpolate(row, col, true, true, false);
        */
/*
            if(pointType == GAP_SURF_GRID)
                setPointType(row, col, GAP_GRID);
            else if(pointType == GAP_NULL || pointType == GAP_NOT_FILLED)
                leastSqFitInterpolate(row, col, true, true, false);
*/
            row += rowInc;
            col += colInc;
        }
    }

    public void fillInsidePoints()
    {
		StsEdgeLoopRadialGridLink link;

        try
        {
            checkInterpolator();

            StsList loopLinks = edgeLoop.getLoopLinks();
            int nLinks = loopLinks.getSize();

            for(int n = 0; n < nLinks; n++)
            {
                link = (StsEdgeLoopRadialGridLink)loopLinks.getElement(n);
                if(link.isConnected(ROW, PLUS)) fillInsidePoints(link, ROW);
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsBlockGrid.fillInsidePoints() failed.",
                e, StsException.WARNING);
        }
    }

    private void fillInsidePoints(StsEdgeLoopRadialGridLink link, int rowOrCol)
    {
        int rowColStart, rowColEnd;
        int nPnts;

        if(link == null) return;
        StsEdgeLoopRadialGridLink nextLink = link.getLink(rowOrCol, PLUS);
        if(nextLink == null) return;
        if(!nextLink.isConnected(rowOrCol, MINUS)) return;  // should be true

        rowColStart = StsMath.ceiling(link.getCrossingRowColF(rowOrCol));
        rowColEnd = StsMath.floor(nextLink.getCrossingRowColF(rowOrCol));
        nPnts = rowColEnd - rowColStart + 1;

        if(nPnts <= 0) return;

        int rowInc, colInc;
        int row, col;

        if(rowOrCol == ROW)
        {
            rowInc = 0;
            colInc = 1;
            row = link.getRowCol(ROW);
            col = rowColStart;
        }
        else
        {
            rowInc = 1;
            colInc = 0;
            row = rowColStart;
            col = link.getRowCol(COL);
        }

        for(int n = 0; n < nPnts; n++)
        {
            byte pointType = getPointType(row, col);
            if(pointType == GAP_NULL || pointType == GAP_NOT_FILLED)
                interpolate(row, col, true, true, false);

            row += rowInc;
            col += colInc;
        }
    }

    StsPolygon getQuadPolygon(int row, int col)
    {
        StsPolygon polygon = new StsPolygon(row, col);
        polygon.type = StsPolygon.SURFACE;
        float[][] cornerPoints = new float[4][3];
        float x = getXCoor(col-colMin);
        float y = getYCoor(row-rowMin);
        cornerPoints[0] = new float[] { x, y, getPointZ(row, col) };
        cornerPoints[1] = new float[] { x, y+yInc, getPointZ(row+1, col) };
        cornerPoints[2] = new float[] { x+xInc, y+yInc, getPointZ(row+1, col+1) };
        cornerPoints[3] = new float[] { x+xInc, y, getPointZ(row, col+1) };
        polygon.setPntNrmls(cornerPoints);
        StsGridRowCol[] gridRowCols = new StsGridRowCol[4];
        gridRowCols[0] = new StsGridRowCol(StsParameters.ROWCOL, row, col);
        gridRowCols[1] = new StsGridRowCol(StsParameters.ROWCOL, row+1, col);
        gridRowCols[2] = new StsGridRowCol(StsParameters.ROWCOL, row+1, col+1);
        gridRowCols[3] = new StsGridRowCol(StsParameters.ROWCOL, row, col+1);
        polygon.gridRowCols = gridRowCols;
        return polygon;
    }
/*
    public void buildNewDomainGeometry(byte zDomain, StsSeismicVelocityModel velocityModel)
    {
        int i = 0, j = 0;
        int row = -1, col = -1;

        try
        {
            float[][] pointsT = getPointsZ();
//            float[][] pointsZ = getAdjPointsZ();
            float xMin = getXMin();
            float xInc = getXInc();
            float yMin = getYMin();
            float yInc = getYInc();

            for(row = pointRowMin; row <= pointRowMax; row++, i++)
            {
                float y = yMin + row*yInc;
                j = 0;
                for(col = pointColMin; col <= pointColMax; col++, j++)
                {
                    float x = xMin + col * xInc;
                    float t = pointsT[i][j];
                    if(t == nullValue)
                        pointsZ[i][j] = nullValue;
                    else
                        pointsZ[i][j] = (float)velocityModel.getZ(x, y, t);
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsEdgeLoop.drawQuadStrips() failed." +
                "row: " + i + " col: " + j, e, StsException.WARNING);
        }
        edgeLoop.buildNewDomainGeometry(zDomain, velocityModel);
    }

    private void buildNewDomainXYSurfaceGeometry(byte zDomain, StsSeismicVelocityModel velocityModel, StsXYSurfaceGridable xyGrid)
    {
        int i = 0, j = 0;
        int row = -1, col = -1;

        try
        {
            float xMin = xyGrid.getXMin();
            float yMin = xyGrid.getYMin();
            float xInc = xyGrid.getXInc();
            float yInc = xyGrid.getYInc();
            int rowMin = xyGrid.getRowMin();
            int colMin = xyGrid.getColMin();
            int rowMax = xyGrid.getRowMax();
            int colMax = xyGrid.getColMax();
            float[][] pointsT = xyGrid.getPointsZ();
//            float[][] pointsZ = xyGrid.getAdjPointsZ();
            for(row = rowMin; row <= rowMax; row++, i++)
            {
                float y = yMin + row*yInc;
                for(col = colMin; col < colMax; col++, j++)
                {
                    float x = xMin + col * xInc;
                    float t = pointsT[i][j];
                    double z = velocityModel.getZ(x, y, t);
                    pointsZ[i][j] = (float)z;
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsEdgeLoop.drawQuadStrips() failed." +
                "row: " + i + " col: " + j, e, StsException.WARNING);
        }
    }
*/

    /** display this surface */
     public void display(StsModel model, StsGLPanel3d glPanel3d, boolean displayFill, boolean displayGrid, String displayMode, boolean displayProperties, StsColor color)
     {
        if(!block.getIsVisible()) return;
        StsPropertyVolume propertyVolume = block.getCurrentPropertyVolume();
         if(displayFill)
             edgeLoop.checkSetPropertyChanged(displayMode);
         if(displayFill || displayGrid)
             if(edgeLoop != null) edgeLoop.display(model, glPanel3d, displayFill, displayGrid, false, false, color, propertyVolume, false);
            // if(edgeLoop != null) edgeLoop.display(model, glPanel3d, displayFill, displayGrid, false, displayProperties, color, propertyVolume, false);
         debugDisplay(currentModel, glPanel3d);
    }

    // called to draw a single gridLine
    public void drawGridLine(int rowOrCol, int rowCol, GL gl)
    {
        if(edgeLoop == null) return;
        boolean displaySurface = surface.getIsVisible() && surface.getDisplayFill();
        StsColor color = displaySurface ? StsColor.BLACK : surface.getStsColor();
  		color.setGLColor(gl);
        edgeLoop.drawGridLines(rowOrCol, rowCol, gl);
    }

    public void displayEdges(StsGLPanel3d glPanel3d)
    {
        StsList edges = getEdges();
        int nEdges = edges.getSize();
        for(int n = 0; n < nEdges; n++)
        {
            StsSurfaceEdge edge = (StsSurfaceEdge)edges.getElement(n);
            edge.display(glPanel3d);
        }
    }

    private void debugDisplay(StsModel model, StsGLPanel3d glPanel3d)
    {
        boolean debug;

        StsBuiltModel builtModel = (StsBuiltModel)model.getCurrentObject(StsBuiltModel.class);
        if(builtModel == null) return;
        StsSurfaceEdge pickedSurfaceEdge = builtModel.pickedSurfaceEdge;
        if(pickedSurfaceEdge == null || pickedSurfaceEdge.getBlockGrid() != this) return;

        if(debugInterpolate) debugDisplayInterpolatePoints(model, glPanel3d);

        debug = model.getBooleanProperty("Block Gap Points");
        if(debug) displayPoints(model, glPanel3d);

        debug = model.getBooleanProperty("debugCellTypes");
        if(debug) displayCellPoints(model, glPanel3d);
    }

    // first point is being interpolated: draw as GAP_FILLED point;
    // draw all subsequent points as GAP_GRID.
    private void debugDisplayInterpolatePoints(StsModel model, StsGLPanel3d glPanel3d)
    {
        StsBlockGridInterpolation.Point point = null;
        GL gl = glPanel3d.getGL();
        int n = -1;

        try
        {
            if(!insideGrid(debugInterpolateRow, debugInterpolateCol)) return;

            float[] xyz;

            gl.glDisable(GL.GL_LIGHTING);
            glPanel3d.setViewShift(gl, 4.0*StsGraphicParameters.gridShift);

            xyz = getXYZorT(debugInterpolateRow, debugInterpolateCol);
            displayPoint(xyz, debugInterpolatePointType, glPanel3d);

            ArrayList points = StsBlockGridInterpolation.getPoints();
            if(points == null) return;

            int nPoints = points.size();
            for(n = 0; n < nPoints; n++)
            {
                point = (StsBlockGridInterpolation.Point)points.get(n);
                xyz = getXYZorT(point.row, point.col);
                byte type = getPointType(point.row, point.col);
                displayPoint(xyz, type, glPanel3d);
            }

            gl.glEnable(GL.GL_LIGHTING);
            glPanel3d.resetViewShift(gl);
        }
        catch(Exception e)
        {
            StsException.outputException("StsBlockGrid.debugDisplayInterpolatePoints() failed." +
			" Interpolation point " + n + ":" + point.toString(),
                e, StsException.WARNING);
        }
    }

    public void displayPoints(StsModel model, StsGLPanel3d glPanel3d)
    {
        try
        {
            float[] point = new float[3];
            GL gl = glPanel3d.getGL();
            gl.glDisable(GL.GL_LIGHTING);

            glPanel3d.setViewShift(gl, 4.0*StsGraphicParameters.gridShift);

            for(int row = pointRowMin; row <= pointRowMax; row++)
            {
                for(int col = pointColMin; col <= pointColMax; col++)
                {
                    point[0] = getXCoor(row, col);
                    point[1] = getYCoor(row, col);
                    byte pointType = getPointType(row, col);
                    point[2] = getPointZ(row, col);
                    displayPoint(point, pointType, glPanel3d);
                }
            }

            gl.glEnable(GL.GL_LIGHTING);
            glPanel3d.resetViewShift(gl);
        }
        catch(Exception e){} // mainDebug, so don't process exception
    }

    // Debug display of point status: BLUE - good grid point (GAP_GRID);
    //                              PURPLE - interpolated point not yet filled (GAP_NOT_FILLED);
    //                               BROWN - interpolated point cant be filled (GAP_CANT_FILLED);
    //                                CYAN - interpolated point (GAP_FILLED);
    //                                 RED - extrapolated point (GAP_EXTRAP);
    //                             MAGENTA - extrapolated point beyond edge: cut (GAP_CUT);
    //                              ORANGE - point is filled hole (GAP_FILL_HOLE).
    //                                GREY - default color.

    static private void displayPoint(float[] point, byte pointType, StsGLPanel3d glPanel3d)
    {
        switch(pointType)
        {
            case GAP_GRID:
                StsGLDraw.drawPoint(point, StsColor.BLUE, glPanel3d, 4, 4, 0.0);
                break;
            case GAP_NOT_FILLED:
                StsGLDraw.drawPoint(point, StsColor.PURPLE, glPanel3d, 4, 4, 0.0);
                break;
            case GAP_CANT_FILL:
                StsGLDraw.drawPoint(point, StsColor.BROWN, glPanel3d, 4, 4, 0.0);
                break;
            case GAP_FILLED:
                StsGLDraw.drawPoint(point, StsColor.CYAN, glPanel3d, 4, 4, 0.0);
                break;
            case GAP_CUT:
                StsGLDraw.drawPoint(point, StsColor.MAGENTA, glPanel3d, 4, 4, 0.0);
                break;
            case GAP_EXTRAP:
                StsGLDraw.drawPoint(point, StsColor.RED, glPanel3d, 4, 4, 0.0);
                break;
            case GAP_FILL_HOLE:
                StsGLDraw.drawPoint(point, StsColor.ORANGE, glPanel3d, 4, 4, 0.0);
                break;
            case GAP_SURF_GRID:
//                StsGLDraw.drawPoint(point, StsColor.WHITE, glPanel3d, 4, 4, 0.0);
                break;
            case GAP_NULL:
                break;
            default:
                StsGLDraw.drawPoint(point, StsColor.GREY, glPanel3d, 4, 4, 0.0);
        }
    }

    public void displayCellPoints(StsModel model, StsGLPanel3d glPanel3d)
    {
        try
        {
            StsColor stsColor;
            float[] point = new float[3];

            if(cellTypeGrid == null)
            {
                StsEdgeLoop edgeLoop = getEdgeLoop();
                if(edgeLoop == null) return;
                cellTypeGrid = edgeLoop.constructCellTypeGrid(this);
                if(cellTypeGrid == null) return;
            }
            GL gl = glPanel3d.getGL();
            gl.glDisable(GL.GL_LIGHTING);

            glPanel3d.setViewShift(gl, 4.0*StsGraphicParameters.gridShift);

            for(int row = pointRowMin; row < pointRowMax; row++)
            {
                for(int col = pointColMin; col < pointColMax; col++)
                {
                    point[0] = getXCoor(row, col);
                    point[1] = getYCoor(row, col);
                    byte cellType = cellTypeGrid.getCellType(row, col);
                    point[2] = getPointZ(row, col);

                    switch(cellType)
                    {
                        case CELL_EMPTY:
                            StsGLDraw.drawPoint(point, StsColor.RED, glPanel3d, 4, 4, 0.0);
                            break;
                        case CELL_FILL:
                            StsGLDraw.drawPoint(point, StsColor.GREEN, glPanel3d, 4, 4, 0.0);
                            break;
                        case CELL_EDGE:
                            StsGLDraw.drawPoint(point, StsColor.YELLOW, glPanel3d, 4, 4, 0.0);
                            break;
                       default:
//                            StsGLDraw.drawPoint(point, StsColor.GREY, glPanel3d, 4, 4, 0.0);
                    }
                }
            }

            gl.glEnable(GL.GL_LIGHTING);
            glPanel3d.resetViewShift(gl);
        }
        catch(Exception e){} // mainDebug, so don't process exception
    }

    public String getLabel()
    {
        return "blockGrid-" + getIndex() + " block: " + block.getLabel() + " surface: " + surface.getName();
    }
    public String toString()
    {
        return "blockGrid-" + getIndex() + " block: " + block.getLabel() + " surface: " + surface.getName();
    }

    public String getPointString(int row, int col)
    {
        byte type = getPointType(row, col);
        if(type == GAP_NULL) return "";
        String typeName = StsParameters.getGapTypeName(type);
        float z = getPointZ(row, col);
        return " " + block.getLabel() + " " + typeName + " " + z;
    }

    public boolean toggleSurfacePickingOn()
    {
        return surface.toggleSurfacePickingOn();
    }
    public void toggleSurfacePickingOff()
    {
        surface.toggleSurfacePickingOff();
    }
    public String getName()
    {
        return surface.getName();
    }
    public StsGridPoint getSurfacePosition(StsMouse mouse, boolean display, StsGLPanel3d glPanel3d)
    {
        return surface.getSurfacePosition(mouse, display, glPanel3d);
    }
    public void setIsVisible(boolean isVisible)
    {
        surface.setIsVisible(isVisible);
    }

    public boolean getIsVisible()
    {
        return surface.getIsVisible();
    }   

    public void propertyTypeChanged()
    {
        edgeLoop.propertyTypeChanged();
    }
/*
    private void debugDrawGridLines(StsWin3d win3d, GL gl)
    {
        StsLinkGridLine gridLine;

        try
        {
            StsSurfaceEdge currentEdge = StsSurfaceEdge.getCurrentEdge();
            if(currentEdge == null || currentEdge.getEdgeGroup() != edgeGroup) return;

            gl.glDisable(GL.GL_LIGHTING);
            win3d.glPanel3d.setViewShift(gl, StsGraphicParameters.gridShift);

            gl.glLineWidth(StsGraphicParameters.gridLineWidth);

            StsSpectrum spectrum = currentModel.getSpectrum("Basic");
            StsColor.setGLColor(gl, spectrum.getColor(edgeGroup.index()%30) );

            for(int row = pointRowMin; row <= pointRowMax; row++)
            {
                gridLine = getFirstGridLine(ROW, row);
                while(gridLine != null)
                {
                    debugDrawGridLine(gridLine, gl);
                    gridLine = gridLine.getNext();
                }
            }

            for(int col = pointColMin; col <= pointColMax; col++)
            {
                gridLine = getFirstGridLine(COL, col);
                while(gridLine != null)
                {
                    debugDrawGridLine(gridLine, gl);
                    gridLine = gridLine.getNext();
                }
            }

            gl.glEnable(GL.GL_LIGHTING);
            win3d.glPanel3d.resetViewShift(gl);
        }
        catch(Exception e){} // mainDebug, so don't process exception
    }

    private void debugDrawGridLine(StsLinkGridLine gridLine, GL gl)
    {
        StsEdgeLoopRadialGridLink endLink;
        float[] xyz;
        int row, col;
        int rowOrCol = NONE;

        try
        {
            rowOrCol = gridLine.getRowOrCol();

            gl.glBegin(GL.GL_LINE_STRIP);

            endLink = gridLine.getFirstLink();
            if(endLink != null)
            {
                xyz = endLink.getXYZ();
                if(xyz != null)  gl.glVertex3fv(xyz);
            }

            int first = gridLine.getFirst();
            int last = gridLine.getLast();

            if(first <= last)
            {
                if(rowOrCol == ROW)
                {
                    row = gridLine.getRowCol();
                    for(col = first; col <= last; col++)
                    {
                        xyz = getXYZ(row, col);
                        if(xyz != null)  gl.glVertex3fv(xyz);
                    }
                }
                else
                {
                    col = gridLine.getRowCol();
                    for(row = first; row <= last; row++)
                    {
                        xyz = getXYZ(row, col);
                        if(xyz != null)  gl.glVertex3fv(xyz);
                    }
                }
            }

            endLink = gridLine.getLastLink();
            if(endLink != null)
            {
                xyz = endLink.getXYZ();
                if(xyz != null)  gl.glVertex3fv(xyz);
            }

            gl.glEnd();
        }
        catch(Exception e)
        {
            int rowCol = gridLine.getRowCol();

            StsException.outputException("StsBlockGrid.debugDrawGridLine() failed. " +
                StsParameters.rowCol(rowOrCol) + rowCol, e, StsException.WARNING);
        }
    }

    public float[][] getLineXYZs(StsLinkGridLine gridLine, float rowColF0, float rowColF1,
                        int lowerBound, int upperBound)
    {
        StsEdgeLoopRadialGridLink firstLink, lastLink;
        float crossingRowColF;

        try
        {
            float firstF = gridLine.getFirstF();
            if(rowColF1 <= firstF) return null;

            float lastF = gridLine.getLastF();
            if(rowColF0 >= lastF) return null;

            int nPnts = upperBound - lowerBound + 1;

            firstLink = gridLine.getFirstLink();
            lastLink = gridLine.getLastLink();

            int first = Math.max(gridLine.getFirst(), lowerBound);
            int last = Math.min(gridLine.getLast(), upperBound);

            float[][] XYZs = new float[nPnts][];

            int rowOrCol = gridLine.getRowOrCol();

            int n = 0;
            if(firstLink != null)
            {
                crossingRowColF = firstLink.getCrossingRowColF(rowOrCol);
                if(crossingRowColF >= lowerBound)
                    XYZs[n++] = firstLink.getXYZ();
            }

            if(rowOrCol == ROW)
            {
                int row = gridLine.getRowCol();
                for(int col = first; col <= last; col++)
                    XYZs[n++] = getXYZ(row, col);
            }
            else if(rowOrCol == COL)
            {
                int col = gridLine.getRowCol();
                for(int row = first; row <= last; row++)
                    XYZs[n++] = getXYZ(row, col);
            }

            if(lastLink != null)
            {
                crossingRowColF = lastLink.getCrossingRowColF(rowOrCol);
                if(crossingRowColF <= upperBound)
                    XYZs[n++] = lastLink.getXYZ();
            }

            if(n != nPnts)
            {
                StsException.systemError("StsBlockGrid.getLineXYZs failed." +
                    " Number of points for gridLine: " + n + " doesn't agree with number expected: " + nPnts);
                return null;
            }
            return XYZs;
        }
        catch(Exception e)
        {
            StsException.outputException("StsBlockGrid.getLineXYZs failed.", e, StsException.WARNING);
            return null;
        }
    }

    public float[][] getLineXYZs(StsLinkGridLine gridLine, float rowColF0, float rowColF1)
    {
        StsEdgeLoopRadialGridLink firstLink, lastLink;

        float firstF = gridLine.getFirstF();
        if(rowColF1 <= firstF) return null;

        float lastF = gridLine.getLastF();
        if(rowColF0 >= lastF) return null;

        int nPnts = 0;

        firstLink = gridLine.getFirstLink();
        if(firstLink != null) nPnts++;

        lastLink = gridLine.getLastLink();
        if(lastLink != null) nPnts++;

        int first = gridLine.getFirst();
        int last = gridLine.getLast();
        nPnts += last-first+1;

        float[][] XYZs = new float[nPnts][];

        int n = 0;
        if(firstLink != null) XYZs[n++] = firstLink.getXYZ();

        int rowOrCol = gridLine.getRowOrCol();

        if(rowOrCol == ROW)
        {
            int row = gridLine.getRowCol();
            for(int col = first; col <= last; col++)
                XYZs[n++] = getXYZ(row, col);
        }
        else if(rowOrCol == COL)
        {
            int col = gridLine.getRowCol();
            for(int row = first; row <= last; row++)
                XYZs[n++] = getXYZ(row, col);
        }

        if(lastLink != null) XYZs[n++] = lastLink.getXYZ();

        return XYZs;
    }

    public float[][] getLineRowColFs(StsLinkGridLine gridLine, float rowColF0, float rowColF1)
    {
        StsEdgeLoopRadialGridLink firstLink, lastLink;
        int row, col, n;

        float firstF = gridLine.getFirstF();
        if(rowColF1 <= firstF) return null;

        float lastF = gridLine.getLastF();
        if(rowColF0 >= lastF) return null;

        int nPnts = 0;

        firstLink = gridLine.getFirstLink();
        if(firstLink != null) nPnts++;

        lastLink = gridLine.getLastLink();
        if(lastLink != null) nPnts++;

        int first = gridLine.getFirst();
        int last = gridLine.getLast();
        nPnts += last-first+1;

        float[][] rowColFs = new float[nPnts][2];

        int rowOrCol = gridLine.getRowOrCol();
        if(rowOrCol == ROW)
        {
            row = gridLine.getRowCol();

            for(n = 0; n < nPnts; n++)
                rowColFs[n][0] = row;

             n = 0;
            if(firstLink != null) rowColFs[n++][1] = firstLink.getCrossingRowColF(ROW);
            for(col = first; col <= last; col++)
                rowColFs[n++][1] = col;
            if(lastLink != null) rowColFs[n++][1] = lastLink.getCrossingRowColF(ROW);
        }
        else if(rowOrCol == COL)
        {
            n = 0;
            col = gridLine.getRowCol();

            for(n = 0; n < nPnts; n++)
                rowColFs[n][1] = col;

            n = 0;
            if(firstLink != null) rowColFs[n++][0] = firstLink.getCrossingRowColF(COL);
            for(row = first; row <= last; row++)
                rowColFs[n++][0] = row;
            if(lastLink != null) rowColFs[n++][0] = lastLink.getCrossingRowColF(COL);
        }

        return rowColFs;
    }
*/
}
