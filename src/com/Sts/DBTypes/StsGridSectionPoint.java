
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System


package com.Sts.DBTypes;

import com.Sts.DB.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

import java.io.*;     

/** This class represents a point on an edge on a section.  The edge could be a StsSurfaceEdge,
 *  an StsWellZone, or an StsGridEdge.  gridRowF and gridColF are float index coordinates for
 *  the global (horizon) grid.  sectionRowF and sectinColF are the float index  coordinates of
 * the section grid.
 */
public class StsGridSectionPoint extends StsObject implements StsCustomSerializable, Cloneable
{
    protected StsPoint point;           /** @param point location of this point                     */
    protected StsSurfaceVertex vertex;  /** @param vertex this point may be on                      */
    protected StsGridRowCol gridRowCol;
    protected StsSectionRowCol[] sectionRowCols;

    transient protected StsEdgeLinkable edge;     /** @param edge edge this point is on  */

    // State flags for connect
    static public final int MINUS = StsParameters.MINUS;
    static public final int PLUS = StsParameters.PLUS;
    static public final int NONE = StsParameters.NONE;
    static public final int PLUS_AND_MINUS = StsParameters.PLUS_AND_MINUS;

    // State flags for rowOrCol
    static final int ROW = StsParameters.ROW;
    static final int COL = StsParameters.COL;
    static final int ROWCOL = StsParameters.ROWCOL;

    static final float nullValue = StsParameters.nullValue;
    static final float largeFloat = StsParameters.largeFloat;

	// this class uses custom serialization and is always explicitly written out, i.e., it
	// is never part of an StsObjectRefList
    public StsGridSectionPoint()
    {
    }

    public StsGridSectionPoint(boolean persistent)
    {
        super(persistent);
    }
    
    public StsGridSectionPoint(StsPoint point, StsSurfaceVertex vertex)
    {
		this(point, vertex.getSectionLine(), vertex);
    }

    public StsGridSectionPoint(StsGridSectionPoint gridSectionPoint, boolean persistent)
    {
        super(persistent);
        StsToolkit.copy(gridSectionPoint, this); // copies only fields from this class and not superclasses, so doesn't include index
    }
/*
    public StsGridSectionPoint(StsPoint point)
    {
        this.point = new StsPoint(point);
        addGridRowCol();
    }
*/
/*
    StsGridSectionPoint(StsPoint wellPoint, StsEdgeLinkable edge, StsWell well, boolean persistent)
    {
        super(persistent);
        this.point = new StsPoint(wellPoint);
        this.edge = edge;
        addGridRowCol(well);
        addSectionRowCols(well);
    }
*/
	// constructor for gridSectionPoint on well at section row
    StsGridSectionPoint(StsPoint point, StsLine line, boolean persistent)
    {
        super(persistent);
        this.point = new StsPoint(point);
        addGridRowCol(line);
        addSectionRowCols(line);
    }

    StsGridSectionPoint(StsPoint wellPoint, StsLine line, StsSurfaceVertex vertex)
    {
        this(wellPoint, line, vertex, true);
    }

    StsGridSectionPoint(StsPoint wellPoint, StsLine line, StsSurfaceVertex vertex, boolean persistent)
    {
        super(persistent);
        this.point = new StsPoint(wellPoint);
        this.vertex = vertex;
		vertex.setSurfacePoint(this);
        addGridRowCol(line);
        addSectionRowCols(line);
    }

    public StsGridSectionPoint(StsPoint point, float rowF, float colF, StsEdgeLinkable edge, StsSurfaceGridable grid, boolean persistent)
    {
        super(persistent);
        this.point = new StsPoint(point);
        this.edge = edge;
        if(grid instanceof StsSection)
        {
            addSectionRowCol(rowF, colF, (StsSection)grid);
            addGridRowCol();
        }
        else
            addGridRowCol(rowF, colF);
    }

    public StsGridSectionPoint(StsGridPoint gridPoint, float rowF, float colF, StsEdgeLinkable edge, StsSurfaceGridable grid, boolean persistent)
    {
        super(persistent);
        this.point = gridPoint.getPoint();
        this.edge = edge;
        if(grid instanceof StsSection)
        {
            addSectionRowCol(rowF, colF, (StsSection)grid);
            addGridRowCol();
        }
        else
            addGridRowCol(rowF, colF);
    }

    public StsGridSectionPoint(StsPoint point, int row, int col, StsEdgeLinkable edge, StsSurfaceGridable grid, boolean persistent)
    {
        super(persistent);
        this.point = point;
        this.edge = edge;
        if(grid instanceof StsSection)
        {
            addSectionRowCol((float)row, (float)col, (StsSection)grid);
            addGridRowCol();
        }
        else
            addGridRowCol((float)row, (float)col);

//rowColFix
//		setRowOrColIndex(grid, ROW, row);
//		setRowOrColIndex(grid, COL, col);
    }

    public StsGridSectionPoint(StsPoint point, StsSurfaceGridable grid, boolean persistent)
    {
        super(persistent);
        this.point = new StsPoint(point);
        if(grid instanceof StsSection) addSectionRowCol((StsSection)grid);
        addGridRowCol();
    }

    public StsGridSectionPoint(StsPoint point)
    {
        this.point = new StsPoint(point);
        addGridRowCol();
    }

    public StsGridSectionPoint(StsPoint point, float rowF, float colF, boolean persistent)
    {
        super(persistent);
        this.point = point;
        addGridRowCol(rowF, colF);
    }

    private StsGridSectionPoint(StsEdgeLinkable edge, boolean persistent)
    {
        super(persistent);
        this.edge = edge;
    }

    public StsGridSectionPoint(StsGridPoint gridPoint, boolean persistent)
    {
        super(persistent);
        addGridRowCol(gridPoint.rowF, gridPoint.colF);
        point = gridPoint.getPoint();
//        float[] xyzd = gridPoint.pxyz;
//        point = StsPoint.constructXYZD(xyzd);

    }

    public StsGridSectionPoint(StsGridPoint gridPoint, StsLine line)
    {
        addGridRowCol(gridPoint.rowF, gridPoint.colF);
        point = gridPoint.getPoint();
        StsSection section = line.getOnSection();
        if(section == null) return;
        addSectionRowCol(section, line);
    }

    public boolean initialize(StsModel model) { return true; }

    public StsGridSectionPoint copy()
    {
        StsPoint newPoint = point.copy();
        return new StsGridSectionPoint(newPoint, gridRowCol.getRowF(), gridRowCol.getColF(), false);
    }

    public static StsGridSectionPoint sectionInterpolate(StsGridSectionPoint sectionPoint0,
                        StsGridSectionPoint sectionPoint1, double ff,
                        StsEdgeLinkable edge, StsSection section, boolean persistent)
    {
        if(ff == 0.0)
            return sectionPoint0;
        else if(ff == 1.0)
            return sectionPoint1;

        StsGridSectionPoint sectionPoint = new StsGridSectionPoint(edge, persistent);
        sectionPoint.sectionInterpolate(sectionPoint0, sectionPoint1, ff, section);
        sectionPoint.addGridRowCol();
        return sectionPoint;
    }

    public static StsGridSectionPoint sectionInterpolate(StsGridSectionPoint sectionPoint0,
                        StsGridSectionPoint sectionPoint1, double ff, StsEdgeLinkable edge, boolean persistent)
    {
        if(ff == 0.0)
            return sectionPoint0;
        else if(ff == 1.0)
            return sectionPoint1;

        StsGridSectionPoint sectionPoint = new StsGridSectionPoint(edge, persistent);
        sectionPoint.sectionInterpolate(sectionPoint0, sectionPoint1, ff);
        sectionPoint.addGridRowCol();
        return sectionPoint;
    }

    public StsGridSectionPoint getClone()
    {
        try
        {
            return (StsGridSectionPoint)this.clone();
        }
        catch(Exception e)
        {
            StsException.outputException("StsGridSectionPoint.getClone() failed.", e, StsException.WARNING);
            return null;
        }
    }

	public void addSectionRowCols()
	{
		if(vertex == null) return;
		addSectionRowCols(vertex.getSectionLine());
	}

    private void addSectionRowCols(StsLine line)
    {
        if(line == null) return;
        StsSection[] connectedSections = line.getAllSections();
        if(connectedSections == null) return;
        int nConnectedSections = connectedSections.length;
        for(int n = 0; n < nConnectedSections; n++)
        {
            StsSection section = connectedSections[n];
            addSectionRowCol(section);
//rowColFix
//			int col = section.getWellIndex(well);
//			setRowOrColIndex(section, COL, col);
        }
    }

	public void setSectionRow()
	{
		if(sectionRowCols == null) return;
		int nSectionRowCols = sectionRowCols.length;
		for(int n = 0; n < nSectionRowCols; n++)
		{
			StsSection section = sectionRowCols[n].getSection();
			if(!section.isVertical()) sectionRowCols[n].setRowOrCol(ROW);
		}
	}

    public StsRowCol getRowCol(StsSurfaceGridable grid)
    {
        if(grid == null)
            return getGridRowCol();
        else if(grid instanceof StsSection)
            return getSectionRowCol((StsSection)grid);
		else if(grid instanceof StsCursorSection)
            //return getGridRowCol();
            return getCursorRowCol((StsCursorSection)grid);
        else
            return getGridRowCol();
    }

    public StsRowCol getRowCol(StsSurfaceGridable grid, byte polygonType)
    {
        if(polygonType == StsPolygon.GRID || polygonType == StsPolygon.SURFACE)
            return getGridRowCol();
        else if(polygonType == StsPolygon.SECTION)
            return getSectionRowCol((StsSection)grid);
		else if(grid instanceof StsCursorSection)
            //return getGridRowCol();
            return getCursorRowCol((StsCursorSection)grid);
        else
            return getGridRowCol();
    }

    public StsGridRowCol getGridRowCol()
    {
        return gridRowCol;
    }

    public StsSectionRowCol getSectionRowCol(StsSection section)
    {
        if(sectionRowCols != null)
        {
            int nRowCols = sectionRowCols.length;
            for(int n = 0; n < nRowCols; n++)
                if(sectionRowCols[n].getSection() == section) return sectionRowCols[n];
        }
        return addSectionRowCol(section);
    }

    /** In most cases, because of textures, rows are vertical and columns are horizontal on cursorSections.
     *  For edgeloops on StsSections and StsCursorSections, rows are horizontal and columns are vertical.
     *  So flip the row col assignment.
     * @param section
     * @return
     */
    public StsRowCol getCursorRowCol(StsCursorSection section)
    {
        float[] xyz = point.getXYZorT();
        float colF = section.getRowCoor(xyz);
        float rowF = section.getColCoor(xyz);
        return new StsCursorSectionRowCol(rowF, colF, section);
    }
    public StsSectionRowCol[] getSectionRowCols() { return sectionRowCols; }

    public void adjustSectionRowCol(StsSection section)
    {
        float guessColF, colF, rowF;

        if(section == null) return;

        StsSectionRowCol sectionRowCol = getSectionRowCol(section);

        int sectionColPosition = getSectionPosition(section);
        if(sectionColPosition == NONE)
        {
            if(sectionRowCol != null)
                guessColF = sectionRowCol.getColF();
            else
                guessColF = -1.0f;

            colF = section.getColF(point, guessColF);
        }
        else if(sectionColPosition == MINUS)
            colF = 0;
        else
            colF = section.getNCols()-1;

        rowF = section.getRowF(point.getZorT());

        if(sectionRowCol == null)
            addSectionRowCol(rowF, colF, section);
        else
        {
            sectionRowCol.setRowF(rowF);
            sectionRowCol.setColF(colF);
        }
    }

    public int getSectionPosition(StsSection section)
    {
        if(vertex == null) return NONE;
        return section.getSidePosition(vertex.getSectionLine());
    }

    public int getRowOrCol(StsSurfaceGridable grid)
    {
        StsRowCol rowCol = getRowCol(grid);
        if(rowCol == null) return NONE;
        return rowCol.getRowOrCol();
    }

    public float getRowF(StsSurfaceGridable grid)
    {
        StsRowCol rowCol = getRowCol(grid);
        if(rowCol == null) return nullValue;
        return rowCol.getRowF();
    }

    public float getColF(StsSurfaceGridable grid)
    {
        StsRowCol rowCol = getRowCol(grid);
        if(rowCol == null) return nullValue;
        return rowCol.getColF();
    }

    public float getGridColF()
    {
        if(gridRowCol == null) return nullValue;
        return gridRowCol.getColF();
    }

    public float getGridRowF()
    {
        if(gridRowCol == null) return nullValue;
        return gridRowCol.getRowF();
    }

    public float[] getXYZorT() { return point.getXYZorT(); }

     public float[] getXYZorT(boolean isDepth) { return point.getXYZorT(isDepth); }

    public StsPoint getPointXYZorT()
    {
        return new StsPoint(getXYZorT());
    }

   public StsPoint getPointXYZorT(boolean isDepth)
    {
        return new StsPoint(getXYZorT(isDepth));
    }

    public StsPoint getPoint() { return point; }
/*
    public StsSection getSection()
    {
        if(sectionLink == null) return null;
        return sectionLink.getSection();
    }
*/
    public float getZ() { return point.getZorT(); }
    public void setEdge(StsEdgeLinkable edge) { this.edge = edge; }
    public StsEdgeLinkable getEdge() { return edge; }
    public StsSurfaceVertex getVertex() { return vertex; }
    public void setVertex(StsSurfaceVertex vertex) { this.vertex = vertex; }

    public int getSide()
    {
        if(edge instanceof StsSurfaceEdge)
            return ((StsSurfaceEdge)edge).getSide();
        else
            return NONE;
    }

    public float getRowColIndexF(StsSurfaceGridable grid, int rowOrCol)
    {
		StsRowCol rowCol = getRowCol(grid);
		if(rowCol == null) return nullValue;
        return rowCol.getRowColIndexF(rowOrCol);
    }

    public int getRowColIndex(StsSurfaceGridable grid, int rowOrCol)
    {
		StsRowCol rowCol = getRowCol(grid);
		if(rowCol == null) return -1;
        return rowCol.getRowCol(rowOrCol);
    }

    public float getCrossingRowColF(StsSurfaceGridable grid)
    {
		StsRowCol rowCol = getRowCol(grid);
		if(rowCol == null) return nullValue;
        return rowCol.getCrossingRowColF();
    }

    public float getCrossingRowColF(StsSurfaceGridable grid, int rowOrCol)
    {
		StsRowCol rowCol = getRowCol(grid);
		if(rowCol == null) return nullValue;
        return gridRowCol.getCrossingRowColF(rowOrCol);
    }


    public void setRowOrColIndex(StsSurfaceGridable grid, int rowOrCol, int index)
    {
		setRowOrColIndex(grid, rowOrCol, index, true);
    }

    public void setRowOrColIndex(StsSurfaceGridable grid, int rowOrCol, int index, boolean checkRedundant)
    {
		StsRowCol rowCol = getRowCol(grid);
		if(rowCol == null) return;
        rowCol.setRowOrColIndex(rowOrCol, index, checkRedundant);
    }

    public void setRowOrColIndexF(StsSurfaceGridable grid, int rowOrCol, float indexF)
    {
		StsRowCol rowCol = getRowCol(grid);
		if(rowCol == null) return;
        rowCol.setRowOrColIndexF(rowOrCol, indexF);
    }
/*
    public void initializeGridInteriorPoint()
    {
        if(gridLink == null) gridLink = new StsGridRowColLink(this, 4, NEED_CONNECT);
        gridRowCol.computeRowColF();
    }
*/
/*
    public void initializeGridEdgePoint()
    {
        if(gridLink == null) gridLink = new StsGridRowColLink(this, 6, DONT_CONNECT);
        gridRowCol.computeRowColF();
    }
*/
/*
    public void computeGridRowColF(int nLinks)
    {
        if(gridLink == null) gridLink = new StsGridRowColLink(this, nLinks, DONT_CONNECT);
        gridRowCol.computeRowColF();
    }
*/
    /** If this point is on a well on this section or on this section,
     *  classInitialize the rowCol info.
     */

    public StsRowCol addRowCol(StsSurfaceGridable grid)
    {
        if(grid instanceof StsSection)
            return addSectionRowCol((StsSection)grid);
        else
            return addGridRowCol();
    }

	public StsSectionRowCol addSectionRowCol(float rowF, float colF, StsSection section)
	{
        StsSectionRowCol sectionRowCol = new StsSectionRowCol(rowF, colF, section);
        addSectionRowCol(sectionRowCol);
        return sectionRowCol;
    }

	public StsRowCol addGridRowCol(float rowF, float colF)
	{
        if(gridRowCol != null)
            StsException.systemError("StsGridSectionPoint.addRowCol() failed: gridRowCol already set.");
        else
            gridRowCol = new StsGridRowCol(rowF, colF);

        return gridRowCol;
    }


	public StsSectionRowCol addSectionRowCol(StsSection section, StsLine line)
	{
		if(section.sectionPatch == null) return null;
        StsSectionRowCol sectionRowCol = new StsSectionRowCol(this, section, line);
        addSectionRowCol(sectionRowCol);
        return sectionRowCol;
    }

	public StsSectionRowCol addSectionRowCol(StsSection section)
	{
		if(section == null || section.sectionPatch == null) return null;
        StsSectionRowCol sectionRowCol = new StsSectionRowCol(this, section);
        addSectionRowCol(sectionRowCol);
        return sectionRowCol;
    }

	public StsRowCol addGridRowCol()
	{
//        if(gridRowCol != null)
//            StsException.systemError("StsGridSectionPoint.addRowCol() failed: gridRowCol already set.");
//        else
        gridRowCol = new StsGridRowCol(point);
        return gridRowCol;
    }

	public StsRowCol addGridRowCol(StsLine line)
	{
//        if(gridRowCol != null)
//            StsException.systemError("StsGridSectionPoint.addRowCol() failed: gridRowCol already set.");
//        else
        gridRowCol = new StsGridRowCol(point);
		StsSection[] sections = line.getAllSections();
		int currentRowOrCol = NONE;
		for(int n = 0; n < sections.length; n++)
		{
			int rowOrCol = sections[n].getRowColAlignment();
			currentRowOrCol = resetRowOrCol(currentRowOrCol, rowOrCol);
		}
		if(currentRowOrCol != NONE) gridRowCol.setRowOrCol(currentRowOrCol);
        return gridRowCol;
    }

	static private int resetRowOrCol(int currentRowOrCol, int newRowOrCol)
	{
		if(currentRowOrCol == ROWCOL || newRowOrCol == ROWCOL)
			return ROWCOL;
		else if(currentRowOrCol == ROW)
		{
			if(newRowOrCol == COL) return ROWCOL;
			else                   return ROW;
		}
		else if(currentRowOrCol == COL)
		{
			if(newRowOrCol == ROW) return ROWCOL;
			else return COL;
		}
		else // currentRowOrCol == NONE
			return newRowOrCol;
	}

    public void addSectionRowCol(StsSectionRowCol sectionRowCol)
    {
        if(sectionRowCols == null)
            sectionRowCols = new StsSectionRowCol[] { sectionRowCol };
        else
        {
            int nSectionRowCols = sectionRowCols.length;
            StsSection section = sectionRowCol.getSection();
            for(int n = 0; n < nSectionRowCols; n++)
            {
                if(sectionRowCols[n].getSection() == section)
				{
				    sectionRowCols[n] = sectionRowCol;
					return;
				}
			/*
                if(sectionRowCols[n].getSection() == section)
                {
                    StsException.systemError("StsGridSectionPoint.addSectionRowCol() failed." +
                        " This section: " + section.getLabel() + " has already been added.");
                    return;
                }
			*/
            }

            StsSectionRowCol[] oldSectionRowCols = sectionRowCols;
            sectionRowCols = new StsSectionRowCol[nSectionRowCols+1];
            System.arraycopy(oldSectionRowCols, 0, sectionRowCols, 0, nSectionRowCols);
            sectionRowCols[nSectionRowCols] = sectionRowCol;
        }
    }

    public void addSectionRowCol(StsSection section, float rowF, float colF, int rowOrCol)
    {
		StsSectionRowCol rowCol = addSectionRowCol(rowF, colF, section);
        rowCol.setRowOrCol(rowOrCol);
    }

    public void addGridRowCol(float rowF, float colF, int rowOrCol)
    {
		StsRowCol rowCol = addGridRowCol(rowF, colF);
        rowCol.setRowOrCol(rowOrCol);
    }


    // compare this point row-col with the next grid row and col away from
    // this comparePoint (same rowOrCol and direction)
    //
    public boolean matchesRowColIndexes(StsSurfaceGridable grid, int[] rowColIndexes)
    {
		StsRowCol rowCol = getRowCol(grid);
        if(rowCol == null) return false;
        return rowCol.matchesRowCol(rowColIndexes);
    }

    public void setPoint(StsPoint point)
    {
        this.point = point;

		if(gridRowCol != null) gridRowCol.recompute(point);
		if(sectionRowCols != null)
		{
			for(int n = 0; n < sectionRowCols.length; n++)
				sectionRowCols[n].recompute(point);
		}
    }

    /** Check if point is on a row or a col */

    public boolean isRowOrCol(StsSurfaceGridable grid, int matchRowOrCol)
    {
		StsRowCol rowCol = getRowCol(grid);
        return rowCol.isRowOrCol(matchRowOrCol);
    }

    public boolean isRowOrCol(StsSurfaceGridable grid)
	{
		StsRowCol rowCol = getRowCol(grid);
		return rowCol.getRowOrCol() != NONE;
	}

	public boolean isRowAndCol(StsSurfaceGridable grid)
	{
		StsRowCol rowCol = getRowCol(grid);
		return rowCol.getRowOrCol() == ROWCOL;
	}

    public  boolean isRow(StsSurfaceGridable grid)
	{
		StsRowCol rowCol = getRowCol(grid);
		return rowCol.isRow();
	}

    public boolean isCol(StsSurfaceGridable grid)
	{
		StsRowCol rowCol = getRowCol(grid);
		return rowCol.isCol();
	}

    public int getRow(StsSurfaceGridable grid)
	{
		StsRowCol rowCol = getRowCol(grid);
		return rowCol.getRow();
	}

    public int getCol(StsSurfaceGridable grid)
	{
		StsRowCol rowCol = getRowCol(grid);
		return rowCol.getCol();
	}

    public int getGridRow()
    {
        return getRow(null);
    }

    public int getGridCol()
    {
        return getCol(null);
    }

    public boolean isSameRowAndCol(StsSurfaceGridable grid, StsGridSectionPoint otherPoint)
    {
		StsRowCol rowCol = getRowCol(grid);
        StsRowCol otherRowCol = otherPoint.getRowCol(grid);
        return rowCol.isSameRowAndCol(otherRowCol);
    }

    public boolean isSameRowAndCol(StsSurfaceGridable grid, int row, int col)
    {
		StsRowCol rowCol = getRowCol(grid);
        return rowCol.isSameRowAndCol(row, col);
    }

    public boolean RowAndCol(StsSurfaceGridable grid)
	{
		StsRowCol rowCol = getRowCol(grid);
		return rowCol.isRowAndCol();
	}

	public int[] getLowerLeftRowCol(StsSurfaceGridable grid)
	{
	    return getLowerLeftRowCol(grid, null);
	}

	public int[] getLowerLeftRowCol(StsSurfaceGridable grid, StsGridSectionPoint otherPoint)
	{
		StsRowCol rowCol, otherRowCol = null;

		rowCol = getRowCol(grid);
		if(rowCol == null) return null;
		if(otherPoint == null) return rowCol.getRowAndCol();
        otherRowCol = otherPoint.getRowCol(grid);
		return rowCol.getLowerLeftRowCol(otherRowCol);
	}

    public void addConnectedEdge(StsEdgeLinkable edge)
    {
        if(vertex == null) vertex = new StsSurfaceVertex(this);
        vertex.addEdge(edge);
    }

    public boolean isFaulted()
    {
        if(edge != null)
        {
            if(!(edge instanceof StsSurfaceEdge)) return false;
            return ((StsSurfaceEdge)edge).isFaulted();
        }

        StsLine line = getSectionLine();
        if(line != null) return line.isFault() && !line.isDyingFault();

        StsException.systemError("StsGridSectionPoint.isFaulted() failed." +
            getLabel() + " is neither on edge or on fault/well.");
        return false;
    }

    public StsLine getSectionLine()
    {
        if(vertex == null) return null;
        return vertex.getSectionLine();
    }

    public boolean isDyingFault()
    {
        StsLine line = getSectionLine();
        if(line == null) return false;
        return line.isDyingFault();
    }

    public boolean isAuxiliary()
    {
        if(edge == null || !(edge instanceof StsSurfaceEdge)) return false;
        return ((StsSurfaceEdge)edge).isAuxiliary();
    }

    public boolean isBoundary()
    {
        if(edge == null || !(edge instanceof StsSurfaceEdge)) return false;
        return ((StsSurfaceEdge)edge).isBoundary();
    }

    /** Return the two edgePoints on either side of this edgePoint relative to the section.
     *  The first boundingEdgePoint has a lesser row/col value relative to this point,
     *  and the second has a greater row/col value. Either or both can be null.
     */
/*
    public StsGridSectionPoint[] getGridBoundingEdgePoints(int rowOrCol, StsHorizonPointLinkedGrid pointLinkedGrid)
    {
        return gridLink.getBoundingEdgePoints(rowOrCol, pointLinkedGrid);
    }
*/
    public StsLine getDyingFault()
    {
        StsLine line = getSectionLine();
        if(line == null) return null;
        if(line.isDyingFault())
            return line;
        else
            return null;
    }

    public boolean onDyingFault()
    {
        StsLine line = getSectionLine();
        if(line == null) return false;
        return line.isDyingFault();
    }

    public boolean onFault()
    {
        StsLine line = getSectionLine();
        if(line == null) return false;
        return line.isFault();
    }

    /** This does a straight-line interpolation. Don't use if the sectionPoints are not on a
     *  planar portion of the section surface.
     */
    private void sectionInterpolate(StsGridSectionPoint sectionPoint0, StsGridSectionPoint sectionPoint1,
                                        double f, StsSection section)
    {
        StsPoint point0, point1;
        try
        {
            point0 = sectionPoint0.getPoint();
            point1 = sectionPoint1.getPoint();
            point = StsPoint.staticInterpolatePoints(point0, point1, (float)f);
            if(section != null) interpolateSectionRowCols(sectionPoint0, sectionPoint1, f);
        }
        catch(Exception e)
        {
            StsException.outputException("StsGridSectionPoint.sectionInterpolate() failed.",
                e, StsException.WARNING);
        }
    }
/*
    private void interpolateSectionRowCols(StsGridSectionPoint point0, StsGridSectionPoint point1,
                                           double f, StsSection section)
    {
        StsSectionRowCol rowCol0 = point0.getSectionRowCol(section);
        StsSectionRowCol rowCol1 = point1.getSectionRowCol(section);
        try
        {
            StsSectionRowCol sectionRowCol = new StsSectionRowCol(rowCol0, rowCol1, f);
            addSectionRowCol(sectionRowCol);
        }
        catch(Exception e)
        {
            StsException.systemError(this, "interpolateSectionRowCols", "failed for section: " + section.getName());
        }
    }
*/
    private void sectionInterpolate(StsGridSectionPoint sectionPoint0, StsGridSectionPoint sectionPoint1, double f)
    {
        StsPoint point0, point1;
        try
        {
            if(point == null) point = new StsPoint(3);

            point0 = sectionPoint0.getPoint();
            point1 = sectionPoint1.getPoint();

            point.interpolatePoints(point0, point1, f);
            interpolateSectionRowCols(sectionPoint0, sectionPoint1, f);
        }
        catch(Exception e)
        {
            StsException.outputException("StsGridSectionPoint.interpolateSectionPoint() failed.",
                e, StsException.WARNING);
        }
    }

    private void interpolateSectionRowCols(StsGridSectionPoint point0, StsGridSectionPoint point1, double f)
    {
        StsSection section = null;

        StsSectionRowCol[] sectionRowCols = point0.getSectionRowCols();
        if(sectionRowCols == null) return;
        try
        {
            for(StsSectionRowCol rowCol0 : sectionRowCols)
            {
                section = rowCol0.getSection();
                StsSectionRowCol rowCol1 = point1.getSectionRowCol(section);
                if(rowCol1 == null) continue;
                StsSectionRowCol sectionRowCol = new StsSectionRowCol(rowCol0, rowCol1, f);
                addSectionRowCol(sectionRowCol);
            }
        }
        catch(Exception e)
        {
            StsException.systemError(this, "interpolateSectionRowCols", "failed for section: " + section.getName());
        }
    }

    private void interpolateGridRowCols(StsGridSectionPoint point0, StsGridSectionPoint point1,
                                           double f)
    {
        StsGridRowCol rowCol0 = point0.getGridRowCol();
        StsGridRowCol rowCol1 = point1.getGridRowCol();
        gridRowCol = new StsGridRowCol(rowCol0, rowCol1, f);
    }

    public float[] getCoordinates()
    {
        return point.v;
    }
/*
    public StsEdgeLoop getEdgeLoop()
    {
        if(edgeLoop != null) return edgeLoop;

        if(edge != null && edge instanceof StsSurfaceEdge)
            return ((StsSurfaceEdge)edge).getEdgeLoop();

        return null;
    }
*/

    public static StsGridSectionPoint[] getTrimmedArray(StsList sectionPointsList)
    {
        int nSectionPoints = sectionPointsList.getSize();
        StsGridSectionPoint[] sectionPoints = new StsGridSectionPoint[nSectionPoints];
        Object[] objects = sectionPointsList.getList();
        for(int n = 0; n < nSectionPoints; n++)
            sectionPoints[n] = (StsGridSectionPoint)objects[n];
        return sectionPoints;
    }


    public void debugDisplayGridPoint(StsGLPanel3d glPanel3d)
    {
//        gridLink.debugDisplayPoint(glPanel3d);
    }

    public static StsList getGridCrossings(StsXYSurfaceGridable grid, StsSection section,
                          StsGridSectionPoint sectionPoint0, StsGridSectionPoint sectionPoint1,
                          StsEdgeLinkable edge, boolean persistent, int endsFlag)
    {
        return getGridCrossings(grid, section, sectionPoint0, sectionPoint1, edge, NONE, persistent, endsFlag);
    }

    public static StsList getGridCrossings(StsGridSectionPoint sectionPoint0, StsGridSectionPoint sectionPoint1,
                                                StsSection section, boolean persistent, int endsFlag)
    {
        return getGridCrossings(null, section, sectionPoint0, sectionPoint1, null, NONE, persistent, endsFlag);
    }

    public static StsList getGridCrossings(StsXYSurfaceGridable grid, StsSection section,
                          StsGridSectionPoint sectionPoint0, StsGridSectionPoint sectionPoint1,
                          StsEdgeLinkable edge, int extrapolate, boolean persistent, int endsFlag)
    {
        StsPoint point0, point1;                /** Input points from sectionPoints */
        float gridRowF0, gridColF0, gridRowF1, gridColF1;               /** float grid coordinates of points */
        StsList gridPoints = null;              /** List of grid crossing points */
        int iStart = 0, iEnd = 0;               /** grid I lines in between the two points */
        int jStart = 0, jEnd = 0;               /** grid J lines in between the two points */
        float dJdI = 0.0f, dIdJ = 0.0f;         /** slope values used in computing intermediate values */
        float gridRowF, gridColF;                           /** float IJ values at the intermediate point */
//        int iConnect = MINUS;                   /** connect directions along J grid line */
//        int  jConnect= MINUS;                   /** connect directions along J grid line */
        float iff = 0.0f, jff = 0.0f;           /** interpolation factors along grid lines */
        StsGridSectionPoint iGridPoint = null;  /** gridPoint along I line */
        StsGridSectionPoint jGridPoint = null;  /** gridPoint along J line */
        StsPoint point;
        boolean iOK, jOK;
        int iInc = 1, jInc = 1;
        int i, j;
        float dI, dJ;

        if(sectionPoint0 == null || sectionPoint1 == null)
            return new StsList(0);

        gridRowF0 = sectionPoint0.getRowF(null);
        gridColF0 = sectionPoint0.getColF(null);

        gridRowF1 = sectionPoint1.getRowF(null);
        gridColF1 = sectionPoint1.getColF(null);

        if(gridRowF1 == gridRowF0 && gridColF1 == gridColF0)
        {
        //    System.out.println("Points are identical.");
            return new StsList(0);
        }

        boolean includeMinusEnd = (endsFlag == MINUS || endsFlag == PLUS_AND_MINUS);
        boolean includePlusEnd = (endsFlag == PLUS || endsFlag == PLUS_AND_MINUS);

        int nEstCrossingPoints = estimateNCrossingPoints(sectionPoint0, sectionPoint1, 1.0f, section);
        gridPoints = new StsList(nEstCrossingPoints, 10); /** List of new points added */

        if(gridRowF0 == gridRowF1)
            iOK = false;
        else
        {
            dJdI = (gridColF1 - gridColF0)/(gridRowF1 - gridRowF0);

            if(gridRowF1 > gridRowF0)
            {
                if(includeMinusEnd)
                    iStart = StsMath.ceiling(gridRowF0);
                else
                    iStart = StsMath.above(gridRowF0);

                if(includePlusEnd)
                    iEnd =  StsMath.floor(gridRowF1);
                else
                    iEnd = StsMath.below(gridRowF1);

                iInc = 1;
                iOK = iEnd >= iStart;
            }
            else // gridRowF1 < gridRowF0
            {
                if(includeMinusEnd)
                    iStart = StsMath.floor(gridRowF0);
                else
                    iStart = StsMath.below(gridRowF0);

                if(includePlusEnd)
                    iEnd = StsMath.ceiling(gridRowF1);
                else
                    iEnd = StsMath.above(gridRowF1);

                iInc = -1;
                iOK = iEnd <= iStart;
            }
        }

        if(gridColF0 == gridColF1)
            jOK = false;
        else
        {
            dIdJ = (gridRowF1 - gridRowF0)/(gridColF1 - gridColF0);

            if(gridColF1 > gridColF0)
            {
                if(includeMinusEnd)
                    jStart = StsMath.ceiling(gridColF0);
                else
                    jStart = StsMath.above(gridColF0);

                if(includePlusEnd)
                    jEnd = StsMath.floor(gridColF1);
                else
                    jEnd = StsMath.below(gridColF1);

                jInc = 1;
                jOK = jEnd >= jStart;
            }
            else // gridColF1 < gridColF0
            {
                if(includeMinusEnd)
                    jStart = StsMath.floor(gridColF0);
                else
                    jStart = StsMath.below(gridColF0);

                if(includePlusEnd)
                    jEnd = StsMath.ceiling(gridColF1);
                else
                    jEnd = StsMath.above(gridColF1);

                jInc = -1;
                jOK = jEnd <= jStart;
            }
        }

        i = iStart;
        j = jStart;

        /** compute location of next point on I grid lines */
        if(iOK)
        {
            gridRowF = (float)i;
            gridColF = gridColF0 + (i - gridRowF0)*dJdI;

            iff = (gridRowF - gridRowF0)/(gridRowF1 - gridRowF0);
            iGridPoint = StsGridSectionPoint.sectionInterpolate(sectionPoint0, sectionPoint1, iff, edge, section, persistent);
            iGridPoint.setRowOrColIndexF(null, ROW, gridRowF);
            if(grid != null) iGridPoint.interpolateBilinearZ(grid);
//            iGridPoint = StsGridSectionPoint.constructor(sectionPoint0, sectionPoint1, iff, section, edge, grid, persistent);
//            iGridPoint.addGridRowCol(gridRowF, gridColF);
        }
        else
            iff = largeFloat;

        /** compute location of next point on J grid lines */
        if(jOK)
        {
            gridRowF =  gridRowF0 + (j - gridColF0)*dIdJ;
            gridColF = (float)j;

            jff = (gridColF - gridColF0)/(gridColF1 - gridColF0);
            jGridPoint = StsGridSectionPoint.sectionInterpolate(sectionPoint0, sectionPoint1, jff, edge, section, persistent);
            jGridPoint.setRowOrColIndexF(null, COL, gridColF);
            if(grid != null) jGridPoint.interpolateBilinearZ(grid);
//            jGridPoint = StsGridSectionPoint.constructor(sectionPoint0, sectionPoint1, jff, section, edge, grid, persistent);
//            jGridPoint.addGridRowCol(gridRowF, gridColF);
        }
        else
            jff = largeFloat;


        while(iOK || jOK)
        {
            if(iff <= jff)
            {
                checkAddGridPoint(gridPoints, iGridPoint);
                i = i + iInc;
                if(iInc > 0 && i <= iEnd || iInc < 0 && i >= iEnd)
                {
                    gridRowF = (float)i;
                    gridColF = (i - gridRowF0)*dJdI + gridColF0;
                    iff = (gridRowF - gridRowF0)/(gridRowF1 - gridRowF0);
                    iGridPoint = StsGridSectionPoint.sectionInterpolate(sectionPoint0, sectionPoint1, iff, edge, section, persistent);
                    iGridPoint.setRowOrColIndexF(null, ROW, gridRowF);
                    if(grid != null) iGridPoint.interpolateBilinearZ(grid);
//                    iGridPoint = StsGridSectionPoint.constructor(sectionPoint0, sectionPoint1, iff, section, edge, grid, persistent);
//                    iGridPoint.addGridRowCol(gridRowF, gridColF);
                }
                else
                {
                    iff = largeFloat;
                    iOK = false;
                }
            }
            else if(jff <= iff)
            {
                checkAddGridPoint(gridPoints, jGridPoint);
                j = j + jInc;
                if(jInc > 0 && j <= jEnd || jInc < 0 && j >= jEnd)
                {
                    gridRowF =  (j - gridColF0)*dIdJ + gridRowF0;
                    gridColF = (float)j;
                    jff = (gridColF - gridColF0)/(gridColF1 - gridColF0);
                    jGridPoint = StsGridSectionPoint.sectionInterpolate(sectionPoint0, sectionPoint1, jff, edge, section, persistent);
                    jGridPoint.setRowOrColIndexF(null, COL, gridColF);
                    if(grid != null) jGridPoint.interpolateBilinearZ(grid);
//                    jGridPoint = StsGridSectionPoint.constructor(sectionPoint0, sectionPoint1, jff, section, edge, grid, persistent);
//                    jGridPoint.addGridRowCol(gridRowF, gridColF);
                }
                else
                {
                    jff = largeFloat;
                    jOK = false;
                }
            }
        }
        return gridPoints;
    }

	public static void adjustPointsToGrid(StsList edgePoints, StsXYSurfaceGridable grid)
	{
		if(edgePoints == null) return;
		int nEdgePoints = edgePoints.getSize();
		for(int n = 0; n < nEdgePoints; n++)
		{
			StsGridSectionPoint edgePoint = (StsGridSectionPoint)edgePoints.getElement(n);
			edgePoint.interpolateBilinearZ(grid);
		}
	}

    public void interpolateBilinearZ(StsXYSurfaceGridable grid)
    {
        float rowF = getGridRowF();
        float colF = getGridColF();
        StsGridPoint gridPoint = new StsGridPoint(point, rowF, colF);
        float z = grid.interpolateBilinearZ(gridPoint, true, true);
        if(z == nullValue)
        {
            StsException.outputException(new StsException(StsException.WARNING,
                "StsGridSectionPoint.interpolateBilinearZ() failed."));
        }
//        point.setT(gridPoint.pxyz[2]);
//        point.setZ(gridPoint.pxyz[3]);

    }

    static public int estimateNCrossingPoints(StsGridSectionPoint point0, StsGridSectionPoint point1, float factor, StsSection section)
    {
        factor = Math.max(1.0f, factor);
        int nEstGridRowCrossings = (int)Math.abs(point1.getRowF(null) - point0.getRowF(null));
        int nEstGridColCrossings = (int)Math.abs(point1.getColF(null) - point0.getColF(null));

        if(section != null)
        {
            int nEstSectionRowCrossings = (int)Math.abs(point1.getRowF(section) - point0.getRowF(section));
            int nEstSectionColCrossings = (int)Math.abs(point1.getColF(section) - point0.getColF(section));
            return (int)(factor*(2 + nEstGridRowCrossings + nEstGridColCrossings + nEstSectionRowCrossings + nEstSectionColCrossings));
        }
        else
            return (int)(factor*(2 + nEstGridRowCrossings + nEstGridColCrossings));
    }

    static private void checkAddGridPoint(StsList gridPoints, StsGridSectionPoint gridPoint)
    {
        StsGridSectionPoint lastPoint = (StsGridSectionPoint)gridPoints.getLast();
        if(lastPoint == null || !lastPoint.isSameRowAndCol(null, gridPoint)) gridPoints.add(gridPoint);
    }

    public void repositionOnSection(StsGridSectionPoint targetSectionPoint)
    {
        point.copyFrom(targetSectionPoint.getPoint());
    }

    static public StsGridSectionPoint[] getPointsVector(StsList pointList, int direction)
    {
        int n, m;

        int nPnts = pointList.getSize();
        StsGridSectionPoint[] pointVector = new StsGridSectionPoint[nPnts];

        if(direction == PLUS)
        {
            for(n = 0; n < nPnts; n++)
                pointVector[n] = (StsGridSectionPoint)pointList.getElement(n);
        }
        else if(direction == MINUS)
        {
            for(n = 0, m = nPnts-1; n < nPnts; n++, m--)
                pointVector[m] = (StsGridSectionPoint)pointList.getElement(n);
        }
        return pointVector;
    }
/*
    public boolean delete()
    {
//        gridLink.delete();
        return true;
    }

    public float gridCompare(int rowOrCol, StsGridSectionPoint otherPoint)
    {
        return gridLink.compare(rowOrCol, otherPoint);
    }

    public int getNeedHasConnectDirection(int rowOrCol)
    {
        return gridLink.getNeedHasConnectDirection(rowOrCol);
    }

    public int getFaultDirection(int rowOrCol)
    {
        return gridLink.getFaultDirection(rowOrCol);
    }

    public int getGridEdgeDirection(int rowOrCol)
    {
        return gridLink.getEdgeDirection(rowOrCol);
    }

    // with this point as center, and the arm towards prevPoint, return the nextPoint
    // obtained by rotating the arm counterclockwise
    public StsGridSectionPoint getGridNextCounterClockwisePoint(StsGridSectionPoint prevPoint)
    {
        return gridLink.getNextCounterClockwisePoint(prevPoint);
    }

    public StsGridSectionPoint getGridNextCounterClockwisePoint(int rowOrCol, int direction)
    {
        return gridLink.getNextCounterClockwisePoint(rowOrCol, direction);
    }

    public StsGridSectionPoint getGridNextCounterClockwisePoint(int linkIndex)
    {
        return gridLink.getNextCounterClockwisePoint(linkIndex);
    }

    // Get next CCW point but restrict to a point in the same grid cell
    public StsGridSectionPoint getGridNextCounterClockwisePolygonPoint(StsGridSectionPoint prevPoint)
    {
        return gridLink.getNextCounterClockwisePolygonPoint(prevPoint);
    }

    public int getGridCounterClockwiseLinkIndex(StsGridSectionPoint prevPoint)
    {
        return gridLink.getCounterClockwiseLinkIndex(prevPoint);
    }

    // with this point as center, and the arm towards prevPoint, return the nextPoint
    public float getGridDistance(int rowOrCol, StsGridSectionPoint nextPoint)
    {
        return gridLink.getDistance(rowOrCol, nextPoint);
    }
*/
    public float getGridDistance(StsSurfaceGridable grid, StsGridSectionPoint nextPoint)
    {
		StsRowCol rowCol, nextRowCol;

		rowCol = getRowCol(grid);
		if(rowCol == null) return nullValue;

		nextRowCol = nextPoint.getRowCol(grid);
		if(nextRowCol == null) return nullValue;

        return rowCol.getGridDistance(nextRowCol);
    }

    public String getRowColLabel(StsSurfaceGridable grid)
    {
	    StsRowCol rowCol = getRowCol(grid);
        return rowCol.toString();
    }

    static public String staticGetLabel(StsGridSectionPoint point)
    {
        if(point == null) return new String("null");
        else return point.getRowColLabel(null);
    }

	public void display(StsGLPanel3d glPanel3d)
	{
		display(glPanel3d, StsColor.WHITE, 3, StsGraphicParameters.vertexShift);
	}

	public void display(StsGLPanel3d glPanel3d, StsColor color, int size, double viewShift)
	{
		if(glPanel3d == null) return;
        if (point == null || point.v == null) return;
        StsGLDraw.drawPoint(point.v, (color==null) ? StsColor.WHITE : color,
                glPanel3d, size, size, viewShift);
    }

    public String toString()
    {
       return toString(point) + " " + toString(gridRowCol);
    }

    static final byte SECTION_EDGE = 1;
	static final byte SURFACE_EDGE = 2;
	static final byte WELL_EDGE = 3;

	public void writeObject(StsDBOutputStream out) throws IllegalAccessException, IOException
	{
		int dbEdgeIndex = 0;
		int dbEdgeType = 0;
		if (edge == null)
		{
			dbEdgeIndex = -1;
			dbEdgeType = -1;
		}
		else if (edge.getClass() == StsSectionEdge.class )
		{
			dbEdgeType = SECTION_EDGE;
		}
		else if (edge.getClass() == StsSurfaceEdge.class )
		{
			dbEdgeType = SURFACE_EDGE;
		}
		else if (edge.getClass() == StsWellZone.class )
		{
			dbEdgeType = WELL_EDGE;
		}
		else
			System.out.println("Unsupported edge in StsGridSectionPoint.write()!");

	    out.writeInt(dbEdgeIndex);
		out.writeInt(dbEdgeType);
	}

	public void readObject(StsDBInputStream in) throws IllegalAccessException, IOException
	{
		int dbEdgeIndex = in.readInt();
		int dbEdgeType = in.readInt();

		edge = null;
		if (dbEdgeIndex != -1)
		{
			switch (dbEdgeType)
			{
				case SECTION_EDGE:
					edge = (StsEdgeLinkable)currentModel.resolveReference(dbEdgeIndex, StsSectionEdge.class);
					break;
				case SURFACE_EDGE:
					edge = (StsEdgeLinkable)currentModel.resolveReference(dbEdgeIndex, StsSurfaceEdge.class);
					break;
				case WELL_EDGE:
					edge = (StsEdgeLinkable)currentModel.resolveReference(dbEdgeIndex, StsWellZone.class);
					break;
				default:
					System.out.println("Unsupported edge in StsGridSectionPoint.read()!");
			}
		}
	}

    public void exportObject(StsDBFileObjectTrader objectTrader) { }
}

