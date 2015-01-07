
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

import java.util.*;

public class StsCellTypeGrid
{
    /** grid this loop is constructed on: an StsSection, or an StsBlockGrid */
    StsSurfaceGridable grid;
    /** Indicates gridType so we know what kind of row and col to get; either SECTION or SURFACE */
    byte gridType;
    /** Edgeloop containing links which bound these cells and partial cells. */
    StsEdgeLoop edgeLoop;
    StsEdgeLoopRadialLinkGrid linkedGrid;


    /** These are cell ranges;  grid ranges extend one more in each direction: rowMax+1, colMax+1 */
    private int rowMin, rowMax, colMin, colMax;
    /** A cell may be EMPTY, FULL, or TRUNCATED */
    private byte[][] cellTypes;

    static final byte CELL_EMPTY = StsParameters.CELL_EMPTY;
    static final byte CELL_FULL = StsParameters.CELL_FULL;
    static final byte CELL_EDGE = StsParameters.CELL_EDGE;
    static final byte CELL_INSIDE_ROW = 3;

    static final String[] CELL_TYPE_STRINGS = new String[] { "CELL_EMPTY", "CELL_FILL", "CELL_EDGE", "CELL_INSIDE_ROW" };
    static public final int MINUS = StsParameters.MINUS;
    static public final int PLUS = StsParameters.PLUS;

    static final int ROW = StsParameters.ROW;
    static final int COL = StsParameters.COL;

    static final int OUTSIDE = 0;
    static final int INSIDE = 2;
    static final int EDGE = 1;

    static final boolean debugCell = false;
    static final int debugCellRow = 62;
    static final int debugCellCol = 48;

    public StsCellTypeGrid(StsSurfaceGridable grid, StsEdgeLoop edgeLoop, StsGridBoundingBox surfaceBoundingBox)
    {
        this.grid = grid;
        this.edgeLoop = edgeLoop;
        if(grid instanceof StsSection)
            gridType = StsPolygon.SECTION;
        else
            gridType = StsPolygon.SURFACE;
        rowMin = StsMath.floor(surfaceBoundingBox.rowMin);
        rowMax = StsMath.ceiling(surfaceBoundingBox.rowMax)-1;
        colMin = StsMath.floor(surfaceBoundingBox.colMin);
        colMax = StsMath.ceiling(surfaceBoundingBox.colMax)-1;
        initializeCellGrid();
        assignLoopCellTypes();
        fillInsideCells();
    }

    public StsCellTypeGrid(StsSurfaceGridable grid, StsEdgeLoop edgeLoop, StsEdgeLoopRadialLinkGrid linkedGrid)
    {
        this.grid = grid;
        this.edgeLoop = edgeLoop;
        this.linkedGrid = linkedGrid;

        // points are from row/col min up to and including max,
        // so cells are from row/col min to max-1
        rowMin = linkedGrid.getRowMin();
        rowMax = linkedGrid.getRowMax();
        colMin = linkedGrid.getColMin();
        colMax = linkedGrid.getColMax();
        initializeCellGrid();
        assignGridCellTypes();
    }

    private void initializeCellGrid()
    {
        if(rowMax < rowMin)
        {
            StsException.systemError("SsCellTypeGrid.consructor() failed. Row range wrong: rowMin " + rowMin + " rowMax " + rowMax);
            return;
        }
        if(colMax < colMin)
        {
             StsException.systemError("SsCellTypeGrid.consructor() failed. Col range wrong: colMin " + colMin + " colMax " + colMax);
             return;
        }

        cellTypes = new byte[rowMax-rowMin+1][colMax-colMin+1];
    }

    private void assignLoopCellTypes()
    {
        StsEdgeLoopRadialGridLink link, nextLink = null;
        Iterator loopIterator = edgeLoop.getLoopGridPointIterator();
        while(loopIterator.hasNext())
        {
            link = (StsEdgeLoopRadialGridLink)loopIterator.next();
            nextLink = link.getNextEdgeLink();
            assignLoopCellType(link, nextLink);
        }
    }

    private void assignLoopCellType(StsEdgeLoopRadialGridLink minusLink, StsEdgeLoopRadialGridLink plusLink)
    {
        StsRowCol minusRowCol = minusLink.point.getRowCol(grid, gridType);
        StsRowCol plusRowCol = plusLink.point.getRowCol(grid, gridType);

        float rowF = minusRowCol.getRowF();
        float colF = minusRowCol.getColF();
        int row = (int)rowF;
        int col = (int)colF;

        float nextRowF = plusRowCol.getRowF();
        float nextColF = plusRowCol.getColF();

        int matchRowOrCol = minusRowCol.getSameRowOrCol(plusRowCol);
        if(matchRowOrCol == ROW)
        {
            if(nextColF > colF)
                row--;
            else if(nextColF < colF && minusRowCol.isCol())
                col--;
            //if(minusRowCol.isCol() && plusRowCol.isCol())
                assignLoopCellType(row, col, CELL_FULL);
            //else
            //    assignLoopCellType(row, col, CELL_EDGE);
        }
        else if(matchRowOrCol == COL)
        {
            if(nextRowF < rowF)
            {
                col--;
                if(minusRowCol.isRow())
                    row--;
            }
            //if(minusRowCol.isRow() && plusRowCol.isRow())
                assignLoopCellType(row, col, CELL_FULL);
            //else
            //    assignLoopCellType(row, col, CELL_EDGE);
        }
        else
        {
            row = Math.min(row, (int)nextRowF);
            col = Math.min(col, (int)nextColF);
            assignLoopCellType(row, col, CELL_EDGE);
        }
    }

    private void assignLoopCellType(int row, int col, byte cellType)
    {
        byte currentCellType = getCellType(row, col);
        if(currentCellType == cellType) return;
        if(currentCellType == CELL_EMPTY)
            setCellType(row, col, cellType);
        else if(currentCellType == CELL_FULL && cellType == CELL_EDGE)
            setCellType(row, col, cellType);
        else if(currentCellType == CELL_EDGE && cellType == CELL_FULL)
            return;
        else
            StsException.systemError(this, "assignLoopCellType", "Undefined operation: assigning " + CELL_TYPE_STRINGS[cellType] + " to " + CELL_TYPE_STRINGS[currentCellType]);
    }

    private void fillInsideCells()
    {
        for(int row = rowMin+1; row < rowMax; row++)
        {
            int firstInsideCol = -1;
            for(int col = colMin+1; col < colMax; col++)
            {
                byte cellType = getCellType(row, col);
                if(cellType == CELL_EMPTY)
                {
                    if(getCellType(row, col-1) != CELL_EMPTY && firstInsideCol == -1)
                        firstInsideCol = col;
                    if(getCellType(row, col+1) != CELL_EMPTY && firstInsideCol != -1)
                    {
                        int lastInsideCol = col;
                        for(int iCol = firstInsideCol; iCol <= lastInsideCol; iCol++)
                            setCellType(row, iCol, CELL_INSIDE_ROW);
                        firstInsideCol = -1;
                    }
                }
            }
        }
        for(int col = colMin+1; col < colMax; col++)
        {
            int firstInsideRow = -1;
            for(int row = rowMin+1; row < rowMax; row++)
            {
                byte cellType = getCellType(row, col);
                if(cellType == CELL_INSIDE_ROW)
                {
                    byte prevRowCellType = getCellType(row-1, col);
                    if(prevRowCellType != CELL_INSIDE_ROW  && prevRowCellType != CELL_EMPTY && firstInsideRow == -1)
                        firstInsideRow = row;
                    byte nextRowCellType = getCellType(row+1, col);
                    if(nextRowCellType != CELL_INSIDE_ROW  && nextRowCellType != CELL_EMPTY && firstInsideRow != -1) 
                    {
                        int lastInsideRow = row;
                        for(int iRow = firstInsideRow; iRow <= lastInsideRow; iRow++)
                            setCellType(iRow, col, CELL_FULL);
                        firstInsideRow = -1;
                    }
                }
            }
        }
    }

    private void assignGridCellTypes()
    {
        StsEdgeLoopRadialGridLink minusLink, plusLink;

        byte[][] pointTypes = new byte[rowMax-rowMin+1][colMax-colMin+1];
        int p00, p01, p10, p11;

        try
        {
            StsList loopLinks = edgeLoop.getLoopLinks();
            int nLinks = loopLinks.getSize();

            for(int n = 0; n < nLinks; n++)
            {
                minusLink = (StsEdgeLoopRadialGridLink)loopLinks.getElement(n);
                if(!minusLink.isRowOrCol()) continue;
                if(minusLink.isRowAndCol())
                {
                    int row = minusLink.getRow();
                    int col = minusLink.getCol();
                    pointTypes[row-rowMin][col-colMin] = INSIDE;
                }
                plusLink = minusLink.getNextInsideConnectedLink(ROW, PLUS);
                while(plusLink != null)
                {
                    assignPointTypes(ROW, minusLink, plusLink, pointTypes);
                    minusLink = plusLink;
                    plusLink = minusLink.getNextInsideConnectedLink(ROW, PLUS);
                }

                plusLink = minusLink.getNextInsideConnectedLink(COL, PLUS);
                while(plusLink != null)
                {
                    assignPointTypes(COL, minusLink, plusLink, pointTypes);
                    minusLink = plusLink;
                    plusLink = minusLink.getNextInsideConnectedLink(COL, PLUS);
                }
            }

            // for a cell to be filled, all four corners must be INSIDE (2), so the sum would be 8
            // for a cell to be empty, all four corners must be outside (0), so the sum would be 0
            // any other summed value means that one of the corners is an edge (just outside), so make the cell an edge cell
            int insideSum;
            for(int row = rowMin; row < rowMax; row++)
            {
                p01 = pointTypes[row-rowMin][0];
                p11 = pointTypes[row+1-rowMin][0];
                for(int col = colMin; col < colMax; col++)
                {
                    p00 = p01;
                    p10 = p11;
                    p01 = pointTypes[row-rowMin][col+1-colMin];
                    p11 = pointTypes[row+1-rowMin][col+1-colMin];
                    insideSum = p00 + p10 + p01 + p11;
                    if(insideSum == 8)
                        setCellType(row, col, CELL_FULL);
                    else if(insideSum > 0)
                        setCellType(row, col, CELL_EDGE);
                    else
                        setCellType(row, col, CELL_EMPTY);
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsEdgeLoop.assignGridCellTypes() failed.",
                e, StsException.WARNING);
        }
        finally
        {
            pointTypes = null;
        }
    }

    private void assignPointTypes(int rowOrCol, StsEdgeLoopRadialGridLink minusLink, StsEdgeLoopRadialGridLink plusLink, byte[][] pointTypes)
    {
		int row, col;
		int colStart, colEnd;
		int rowStart, rowEnd;

		try
		{
			if(rowOrCol == ROW)
			{
				row = minusLink.getRow();
                colStart = StsMath.ceiling(minusLink.getColF());
				colEnd = StsMath.floor(plusLink.getColF());
                if(!minusLink.isCol())
                    pointTypes[row-rowMin][colStart-1-colMin] = EDGE;
				for(col = colStart; col <= colEnd; col++)
					pointTypes[row-rowMin][col-colMin] = INSIDE;
                if(!plusLink.isCol())
                    pointTypes[row-rowMin][colEnd+1-colMin] = EDGE;
            }
			else if(rowOrCol == COL)
			{
				col = minusLink.getCol();
				rowStart = StsMath.ceiling(minusLink.getRowF());
				rowEnd = StsMath.floor(plusLink.getRowF());
                if(!minusLink.isRow())
                    pointTypes[rowStart-1-rowMin][col-colMin] = EDGE;
                for(row = rowStart; row <= rowEnd; row++)
					pointTypes[row-rowMin][col-colMin] = INSIDE;
                if(!plusLink.isRow())
                    pointTypes[rowEnd+1-rowMin][col-colMin] = EDGE;
            }
		}
        catch(Exception e)
        {
            StsException.outputException("StsEdgeLoop.assignPointTypes() failed.",
                e, StsException.WARNING);
        }
    }

    private void setCellType(int row, int col, byte cellType)
    {
        if(!isInsideCellGrid(row, col)) return;

        if(debugCell && row == debugCellRow && col == debugCellCol)
        {
            StsException.systemDebug(" StsCellTypeGrid.setCellType() cellType[" +
                row + "][" + col + "] set to: " + StsParameters.cellTypeString(cellType));
        }
        cellTypes[row-rowMin][col-colMin] = cellType;
    }

    public int getRowMin() { return rowMin; }
    public int getRowMax() { return rowMax; }
    public int getColMin() { return colMin; }
    public int getColMax() { return colMax; }

    public byte getCellType(int[] rowCol)
    {
        return getCellType(rowCol[0], rowCol[1]);
    }

    public byte getCellType(int row, int col)
    {
        if(!isInsideCellGrid(row, col)) return CELL_EMPTY;
        return cellTypes[row-rowMin][col-colMin];
    }

    public boolean isInsideCellGrid(int row, int col)
    {
        return row >= rowMin && row <= rowMax && col >= colMin && col <= colMax;
    }

    public byte[][] getCellTypes()
    {
        return cellTypes;
    }
}
