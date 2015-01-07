
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
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.util.*;

public class StsLinkGrid
{
    protected int rowMin;
    protected int rowMax;
    protected int colMin;
    protected int colMax;
    protected StsEdgeLoopRadialGridLink[] rowLinks;
    protected StsEdgeLoopRadialGridLink[] colLinks;
    protected StsSurfaceGridable grid;
    protected byte gridType; // SECTION or HORIZON

    static final int ROW = StsParameters.ROW;
    static final int COL = StsParameters.COL;
    static final int ROWCOL = StsParameters.ROWCOL;
    static final int NONE = StsParameters.NONE;

    static final int PLUS = StsParameters.PLUS;
    static final int MINUS = StsParameters.MINUS;

    static public final int ROW_PLUS = StsParameters.ROW_PLUS;
    static public final int COL_PLUS = StsParameters.COL_PLUS;
    static public final int ROW_MINUS = StsParameters.ROW_MINUS;
    static public final int COL_MINUS = StsParameters.COL_MINUS;

    static public final byte SECTION = StsParameters.SECTION;
    static public final byte HORIZON = StsParameters.HORIZON;

    static final int nullInteger = StsParameters.nullInteger;
    static final float largeFloat = StsParameters.largeFloat;
    static final float roundOff = StsParameters.roundOff;

    public StsLinkGrid()
    {
    }

    public StsLinkGrid(StsSurfaceGridable grid)
    {
        this.grid = grid;

        if(grid instanceof StsSection)
            gridType = SECTION;
        else
            gridType = HORIZON;

        rowMin = 0;
        rowMax = grid.getNRows()-1;
        colMin = 0;
        colMax = grid.getNCols()-1;

        rowLinks = new StsEdgeLoopRadialGridLink[rowMax-rowMin+1];
        colLinks = new StsEdgeLoopRadialGridLink[colMax-colMin+1];
    }

    public StsLinkGrid(StsGridBoundingBox boundingBox)
    {
        rowMin = boundingBox.rowMin;
        rowMax = boundingBox.rowMax;
        colMin = boundingBox.colMin;
        colMax = boundingBox.colMax;

        int nRows = rowMax-rowMin+1;
        if(nRows < 0 || nRows > 10000) return;

        int nCols = colMax-colMin+1;
        if(nCols < 0 || nCols > 10000) return;

        rowLinks = new StsEdgeLoopRadialGridLink[nRows];
        colLinks = new StsEdgeLoopRadialGridLink[nCols];
    }

    public StsLinkGrid(StsSurfaceGridable grid, StsList edgeLinks)
    {
		this(grid);

        int nLinks = edgeLinks.getSize();
        for(int n = 0; n < nLinks; n++)
        {
            StsEdgeLoopRadialGridLink link = (StsEdgeLoopRadialGridLink)edgeLinks.getElement(n);
            addLink(link);
        }
    }

// Accessors

	public StsSurfaceGridable getGrid() { return grid; }
    public int getRowMin() { return rowMin; }
    public int getRowMax() { return rowMax; }
    public int getColMin() { return colMin; }
    public int getColMax() { return colMax; }

    public void addLink(StsEdgeLoopRadialGridLink link)
    {
        addLink(link, true);
    }

    public void addLink(StsEdgeLoopRadialGridLink link, boolean checkConnect)
    {
        int rowOrCol = link.getRowOrCol();
        if(rowOrCol == NONE) return;
        if(rowOrCol == ROW || rowOrCol == ROWCOL)
            insertLink(ROW, link, checkConnect);
        if(rowOrCol == COL || rowOrCol == ROWCOL)
            insertLink(COL, link, checkConnect);
    }
/*
    public void addPoint(StsGridPointLinkable gridPoint)
    {
        StsPoint point = gridPoint.getPoint();
        float rowF = getRowCoor(point);
        float colF = getColCoor(point);
        StsEdgeLoopRadialGridLink link = new StsEdgeLoopRadialGridLink(gridPoint, rowF, colF);
        addLink(link);
    }
*/
/*
    public float getRowCoor(StsGridSectionPoint gridPoint)
    {
        return grid.getRowCoor(gridPoint.getPoint().v);
    }

    public float getColCoor(StsGridSectionPoint gridPoint)
    {
        return grid.getColCoor(gridPoint.getPoint().v);
    }
*/
    public int getGridMin(int rowOrCol)
    {
        return rowOrCol == ROW ? colMin : rowMin;
    }

    public int getGridMax(int rowOrCol)
    {
        return rowOrCol == ROW ? colMax : rowMax;
    }

    public boolean hasRow(int row)
    {
        if(rowLinks == null) return false;
        if(row < rowMin || row > rowMax) return false;
        return rowLinks[row - rowMin] != null;
    }

    public boolean hasCol(int col)
    {
        if(colLinks == null) return false;
        if(col < colMin || col > colMax) return false;
        return colLinks[col - colMin] != null;
    }

    public StsEdgeLoopRadialGridLink getFirstRowLink(int row)
    {
        if(rowLinks == null) return null;
        if(row < rowMin || row > rowMax) return null;
        return rowLinks[row - rowMin];
    }

    public StsEdgeLoopRadialGridLink getFirstColLink(int col)
    {
        if(colLinks == null) return null;
        if(col < colMin || col > colMax) return null;
        return colLinks[col - colMin];
    }

    public StsEdgeLoopRadialGridLink getFirstRowColLink(int rowOrCol, int rowCol)
    {
        try
        {
            if(rowOrCol == ROW)
                return rowLinks[rowCol - rowMin];
            else if(rowOrCol == COL)
                return colLinks[rowCol - colMin];
            else
                return null;
        }
        catch(Exception e)
        {
            StsException.outputException("StsLinkGrid.getFirstRowColLink() failed." + "\n" +
				"Grid: " + grid.getLabel(),
                e, StsException.WARNING);
            return null;
        }
    }

    public StsEdgeLoopRadialGridLink getFirstRowColLink(int rowOrCol, StsEdgeLoopRadialGridLink link)
    {
        if(rowOrCol == NONE) return null;
        int rowCol = link.getRowCol(rowOrCol);
        return getFirstRowColLink(rowOrCol, rowCol);
    }

     public void setFirstRowColLink(int rowOrCol, int rowCol, StsEdgeLoopRadialGridLink link)
     {
        if(rowOrCol == ROW)
            rowLinks[rowCol-rowMin] = link;
        else if(rowOrCol == COL)
            colLinks[rowCol-colMin] = link;
     }

     public void setFirstRowColLink(int rowOrCol, StsEdgeLoopRadialGridLink link)
     {
        try
        {
            if(rowOrCol == ROW)
            {
//rowColFix                int row = StsMath.roundOffInteger(link.getRowF());
                rowLinks[link.getRow()-rowMin] = link;
            }
            else if(rowOrCol == COL)
            {
//rowColFix                int col = StsMath.roundOffInteger(link.getColF());
                colLinks[link.getCol()-colMin] = link;
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsLinkGrid.setFirstRowColLink() failed." + "\n" +
				"Grid: " + grid.getLabel(),
                e, StsException.WARNING);
        }
    }

    public void setStartLink(StsEdgeLoopRadialGridLink link)
    {
        int rowOrCol = link.getRowOrCol();
        if(rowOrCol == ROW)
        {
            int row = (int)link.getRowF();
            rowLinks[row-rowMin] = link;
            return;
        }
        else if(rowOrCol == COL)
        {
            int col = (int)link.getColF();
            colLinks[col-colMin] = link;
            return;
        }
        else
        {
            StsException.systemError("StsLinkLinkedGrid.setStartLink() failed. " +
                link.toString() + " is not row or column: " + StsParameters.rowCol(rowOrCol) + "\n" +
				"Grid: " + grid.getLabel());
            return;
        }
    }
/*
	public void insertPoint(StsGridSectionPoint point)
	{
		StsEdgeLoopRadialGridLink link = new StsEdgeLoopRadialGridLink(grid, point);
		insertLink(link);
	}
*/
    private void insertLink(StsEdgeLoopRadialGridLink link)
    {
        if(link.isRowOrCol(ROW))
            insertLink(ROW, link, true);
        if(link.isRowOrCol(COL))
            insertLink(COL, link, true);
    }

    public void orderLinks()
    {
        for(int row = rowMin; row <= rowMax; row++)
            if(hasRow(row)) orderLinks(ROW, row);

        for(int col = colMin; col <= colMax; col++)
            if(hasCol(col)) orderLinks(COL, col);
    }

    private void orderLinks(int rowOrCol, int rowCol)
    {
        StsEdgeLoopRadialGridLink link;

        ArrayList minusLinks = new ArrayList();
        ArrayList plusLinks = new ArrayList();
        ArrayList midLinks = new ArrayList();

        Enumeration enumeration = getEnum(rowOrCol, rowCol);
        while(enumeration.hasMoreElements())
        {
            link = (StsEdgeLoopRadialGridLink)enumeration.nextElement();
            int direction = link.getInsideDirection(rowOrCol);
            if(direction == MINUS)
                minusLinks.add(link);
            else if(direction == PLUS)
                plusLinks.add(link);
            else
                midLinks.add(link);
        }

        Comparator comparator = new LinkComparator(rowOrCol);
        Collections.sort(minusLinks, comparator);
        Collections.sort(plusLinks, comparator);
        Collections.sort(midLinks, comparator);

        int nMinusLinks = minusLinks.size();
        int nPlusLinks = plusLinks.size();
        int nMidLinks = midLinks.size();
/*
        if(nMinusLinks != nPlusLinks)
        {
            StsException.systemError("StsLinkGrid.orderLinks() failed." +
                " Number of minus and plus links don't compare for: " +
                StsParameters.rowCol(rowOrCol) + " " + rowCol + "\n" +
				"Grid: " + grid.getLabel());
            return;
        }
*/
        StsEdgeLoopRadialGridLink minusLink, plusLink;
        link = null;

        for(int p = 0; p < nPlusLinks; p++)
        {
            plusLink = (StsEdgeLoopRadialGridLink)plusLinks.get(p);

            // get first minusLink which matches this edgeLoop
            minusLink = null;
            Iterator minusLinkIter = minusLinks.iterator();
            while(minusLinkIter.hasNext())
            {
                StsEdgeLoopRadialGridLink mLink = (StsEdgeLoopRadialGridLink)minusLinkIter.next();
                if(mLink.isSameEdgeLoop(plusLink))
                {
                    minusLink = mLink;
                    minusLinks.remove(mLink);
                    break;
                }
            }

            if(minusLink == null)
            {
                StsException.systemError("StsLinkGrid.orderLinks() failed." +
                    " Couldn't find minus link matching plus link edgeLoop: " + plusLink.toString() + "\n" +
				    "Grid: " + grid.getLabel());
                continue;
            }

            if(link == null)
            {
                setFirstRowColLink(rowOrCol, rowCol, plusLink);
                plusLink.setLink(rowOrCol, MINUS, null);
            }
            else
                link.connectToNextLink(rowOrCol, plusLink, grid);

            link = plusLink;

            StsEdgeLoopRadialGridLink[] midLinksArray =  (StsEdgeLoopRadialGridLink[])midLinks.toArray(new StsEdgeLoopRadialGridLink[0]);
            nMidLinks = midLinksArray.length;

            for(int n = 0; n < nMidLinks; n++)
            {
                StsEdgeLoopRadialGridLink midLink = midLinksArray[n];
                if(isBetween(rowOrCol, midLink, plusLink, minusLink))
                {
                    link.connectToNextLink(rowOrCol, midLink, grid);
                    link = midLink;
                    midLinks.remove(midLink);
                }
            }
            link.connectToNextLink(rowOrCol, minusLink, grid);
            link = minusLink;
        }
        if(link != null) link.setLink(rowOrCol, PLUS, null);
    }
/*
    public void simpleOrderLinks()
    {
        for(int row = rowMin; row <= rowMax; row++)
            if(hasRow(row)) simpleOrderLinks(ROW, row);

        for(int col = colMin; col <= colMax; col++)
            if(hasCol(col)) simpleOrderLinks(COL, col);
    }

    private void simpleOrderLinks(int rowOrCol, int rowCol)
    {
        ArrayList links = new ArrayList();

        Enumeration enum = getEnum(rowOrCol, rowCol);
        while(enum.hasMoreElements())
            links.add((StsEdgeLoopRadialGridLink)enum.nextElement());

        Comparator comparator = new LinkComparator(rowOrCol);
        Collections.sort(links, comparator);

        StsEdgeLoopRadialGridLink link, lastLink = null;
        int nLinks = links.size();
        int plusLinkIndex = StsEdgeLoopRadialGridLink.getLinkIndex(rowOrCol, PLUS);
        int minusLinkIndex = StsEdgeLoopRadialGridLink.getLinkIndex(rowOrCol, MINUS);
        for(int n = 0; n < nLinks; n++)
        {
            link = (StsEdgeLoopRadialGridLink)links.get(n);
            if(lastLink == null)
            {
                setFirstRowColLink(rowOrCol, rowCol, link);
                link.setLink(minusLinkIndex, null);
            }
            else
            {
                lastLink.setLink(plusLinkIndex, link);
                link.setLink(minusLinkIndex, lastLink);
                lastLink.setConnectType(plusLinkIndex, StsEdgeLoopRadialGridLink.FAR_CONNECTED);
                link.setConnectType(minusLinkIndex, StsEdgeLoopRadialGridLink.FAR_CONNECTED);
            }
            lastLink = link;
        }
        if(lastLink != null) lastLink.setLink(plusLinkIndex, null);
    }
*/
    // link.crossingRowCol must be between prevLink and nextLink and
    // all 3 links must belong to the same edgeLoop
    private boolean isBetween(int rowOrCol, StsEdgeLoopRadialGridLink link, StsEdgeLoopRadialGridLink prevLink, StsEdgeLoopRadialGridLink nextLink)
    {
        if(!link.isSameEdgeLoop(prevLink)) return false;
        if(!link.isSameEdgeLoop(nextLink)) return false;
        float rowColF = link.getCrossingRowColF(rowOrCol);
        return isBetween(rowOrCol, rowColF, prevLink, nextLink);
    }

    // link.crossingRowCol must be between prevLink and nextLink on row or col
    private boolean isBetween(int rowOrCol, float rowColF, StsEdgeLoopRadialGridLink prevLink, StsEdgeLoopRadialGridLink nextLink)
    {
        if(rowOrCol != ROW && rowOrCol != COL)
        {
            StsException.systemError("StsEdgeLoopRadialGridLink.isBetween() failed." +
                " rowOrCol argument is " + StsParameters.rowCol(rowOrCol) +
                ". It must be ROW or COL." + "\n" +
				"Grid: " + grid.getLabel());
            return false;
        }

        float prevRowColF = prevLink.getCrossingRowColF(rowOrCol);
        if(rowColF < prevRowColF) return false;
        float nextRowColF = nextLink.getCrossingRowColF(rowOrCol);
        if(rowColF > nextRowColF) return false;

        return true;
    }
   /** LinkComparator compares two links along a row or col sorting them by:
     *  1) loop they are on;
     *  2) increasing order of crossing col or row.
     *  3) z values.
     */
    public final class LinkComparator implements Comparator
    {
        int rowOrCol;

        LinkComparator(int rowOrCol)
        {
            this.rowOrCol = rowOrCol;
        }

        public int compare(Object o1, Object o2)
        {
            StsEdgeLoopRadialGridLink p1 = (StsEdgeLoopRadialGridLink)o1;
            StsEdgeLoopRadialGridLink p2 = (StsEdgeLoopRadialGridLink)o2;
            // sort first by loop; loops are indexed arbitrarily; links will be sorted into this same sequence
            StsEdgeLoop loop1 = p1.getEdgeLoop();
            StsEdgeLoop loop2 = p2.getEdgeLoop();
            if(loop1 != null && loop2 != null)
            {
                int index1 = loop1.getIndex();
                int index2 = loop2.getIndex();
                if(index1 > index2) return 1;
                if(index1 < index2) return -1;
            }
            // if on same loop, sort by crossing col if on row, or by crossing row if on col
            float crossingRowCol1 = p1.getCrossingRowColF(rowOrCol);
            float crossingRowCol2 = p2.getCrossingRowColF(rowOrCol);

            if(crossingRowCol1 > crossingRowCol2) return 1;
            else if(crossingRowCol1 < crossingRowCol2) return -1;
            // if at same row and col, sort by z
            float z1 = p1.getPoint().getZ();
            float z2 = p2.getPoint().getZ();
            if(z1 > z2) return 1;
            if(z1 < z2) return -1;

            StsException.systemError(this, "compare", " points are identical: " + p1.toString() + " " + p2.toString() + "\n" + "Grid: " + grid.getLabel());
            return 0;
        }
    }

    // temporarily insert newLink at beginning of row or col list
    private void insertLink(int rowOrCol, StsEdgeLoopRadialGridLink newLink, boolean checkConnect)
    {
        if(rowOrCol == NONE) return;
        if(checkConnect && !newLink.needsConnect(rowOrCol)) return;
        int rowCol = newLink.getRowCol(rowOrCol);

        StsEdgeLoopRadialGridLink link = getFirstRowColLink(rowOrCol, rowCol);
        setFirstRowColLink(rowOrCol, newLink);
        if(link == null) return;
        newLink.initialSetLink(rowOrCol, PLUS, link);
    }

    /** Add a link in the corresponding ROW and or COL.  Reorder when all links added.
     *  The node is at the intersection of row and col grid lines (BOTH), or on
     *  a row between col lines (ROW) or on a col between row lines (COL).
     *  If BOTH, it is inserted into both row and col link lists.  If ROW or COL,
     *  it is inserted only into the corresponding list.
     */
    public void insertOrderedLink(StsEdgeLoopRadialGridLink link)
    {
        int rowOrCol = link.getRowOrCol();
        if(rowOrCol == ROW || rowOrCol == ROWCOL)
            insertOrderedLink(ROW, link);
        if(rowOrCol == COL || rowOrCol == ROWCOL)
            insertOrderedLink(COL, link);
    }

    public StsEdgeLoopRadialGridLink getRowLink(int row)
    {
        if(rowLinks == null) return null;
        if(row >= rowMin && row <= rowMax)
            return rowLinks[row - rowMin];
        else
            return null;
    }

    public StsEdgeLoopRadialGridLink getColLink(int col)
    {
        if(colLinks == null) return null;
        if(col >= colMin && col <= colMax)
            return colLinks[col - colMin];
        else
            return null;
    }

    public int getColMin(int row)
    {
        if(row < rowMin || row > rowMax) return nullInteger;

        StsEdgeLoopRadialGridLink link = rowLinks[row-rowMin];
        if(link == null) return nullInteger;
        return StsMath.ceiling(link.getColF());
    }

    public int getColMax(int row)
    {
        if(row < rowMin || row > rowMax) return nullInteger;

        StsEdgeLoopRadialGridLink link = rowLinks[row-rowMin];
        if(link == null) return nullInteger;

        int linkIndex = StsEdgeLoopRadialGridLink.ROW_PLUS;
        StsEdgeLoopRadialGridLink nextLink;
        while((nextLink = link.getLink(linkIndex)) != null)
            link = nextLink;

        return StsMath.floor(link.getColF());
    }

    public int getRowMin(int col)
    {
        if(col < colMin || col > colMax) return nullInteger;

        StsEdgeLoopRadialGridLink link = colLinks[col-colMin];
        if(link == null) return nullInteger;
        return StsMath.ceiling(link.getRowF());
    }

    public int getRowMax(int col)
    {
        if(col < colMin || col > colMax) return nullInteger;

        StsEdgeLoopRadialGridLink link = colLinks[col-colMin];
        if(link == null) return nullInteger;

        int linkIndex = StsEdgeLoopRadialGridLink.COL_PLUS;
        StsEdgeLoopRadialGridLink nextLink;
        while((nextLink = link.getLink(linkIndex)) != null)
            link = nextLink;

        return StsMath.floor(link.getRowF());
    }


    public boolean delete(int rowOrCol, int rowCol, StsEdgeLoopRadialGridLink deleteLink)
    {
        StsEdgeLoopRadialGridLink firstLink, nextLink;

        int row = deleteLink.getRowCol(ROW);
        if(row >= 0)
        {
            firstLink = getFirstRowColLink(ROW, row);
            if(firstLink == deleteLink)
            {
                nextLink = deleteLink.getLink(ROW, PLUS);
                setFirstRowColLink(ROW, row, nextLink);
            }
        }
        int col = deleteLink.getRowCol(COL);
        if(col >= 0)
        {
            firstLink = getFirstRowColLink(COL, col);
            if(firstLink == deleteLink)
            {
                nextLink = deleteLink.getLink(COL, PLUS);
                setFirstRowColLink(COL, col, nextLink);
            }
        }

        deleteLink.delete();
        return true;
    }

    public boolean isLinkAt(int row, int col)
    {
        Enumeration enumeration = getEnum(ROW, row);

        while(enumeration.hasMoreElements())
        {
            StsEdgeLoopRadialGridLink link = (StsEdgeLoopRadialGridLink)enumeration.nextElement();
            float colF = link.getColF();
            if(colF == (float)col) return true;
        }
        return false;
    }

    // use this method only when adding links on grid locations
    public void insertOrderedLink(int rowOrCol, StsEdgeLoopRadialGridLink newLink)
    {
        float newRowColF, prevRowColF, nextRowColF;
        int linkIndex;
        StsEdgeLoopRadialGridLink prevLink = null, nextLink = null;

        try
        {
            if(newLink == null) return;
            if(newLink.isLinked(rowOrCol, PLUS) && newLink.isLinked(rowOrCol, MINUS)) return;

            if(rowOrCol == NONE) return;
            int rowCol = newLink.getRowCol(rowOrCol);

            nextLink = getFirstRowColLink(rowOrCol, rowCol);
            if(nextLink == null)
            {
                setFirstRowColLink(rowOrCol, newLink);
                return;
            }

            linkIndex = StsEdgeLoopRadialGridLink.getLinkDirection(rowOrCol, PLUS);

            prevLink = nextLink;
            nextLink = prevLink.getLink(linkIndex);

            while(nextLink != null)
            {
                if(isBetween(rowOrCol, newLink, prevLink, nextLink))
                {
                    newLink.connectToLinks(rowOrCol, prevLink, nextLink, grid);
//                    prevLink.setRowColConnectedLink(rowOrCol, newLink);
//                    newLink.setRowColConnectedLink(rowOrCol, nextLink);
                    return;
                }
                prevLink = nextLink;
                nextLink = prevLink.getLink(linkIndex);
            }
            StsException.systemError("StsLinkGrid.insertOrderedLink() failed. " +
                "Could not insert link: " + newLink.toString() + "\n" +
				"Grid: " + grid.getLabel());
        }
        catch(Exception e)
        {
            StsException.outputException("StsLinkGrid.insertOrderedLink() failed. " +
                "Could not insert link: " + newLink.toString() + "\n" +
				"Grid: " + grid.getLabel(), e, StsException.WARNING);
        }
    }

    /** Returns an enumeration of the components of this StsList.
     * @return  an enumeration of the components of this list.
     */

    public final Enumeration getEnum(int rowOrCol, int rowCol)
    {
        if(rowOrCol == ROW)
            return rowEnum(rowCol);
        else
            return colEnum(rowCol);
    }

    public final Enumeration rowEnum(int row)
    {
        if(row < rowMin || row > rowMax)
            return new LineEnumerator();
        else
	        return new LineEnumerator(ROW, row);
    }

    public final Enumeration colEnum(int col)
    {
        if(col < colMin || col > colMax)
            return new LineEnumerator();
        else
	        return new LineEnumerator(COL, col);
    }

    /** StsListEnumerator iterates over entries in this row */
    final class LineEnumerator implements Enumeration
    {
        int index = -1;
        StsEdgeLoopRadialGridLink nextLink = null;
        int linkIndex;

        LineEnumerator()
        {
        }

        LineEnumerator(int rowOrCol, int index)
        {
            this.index = index;

            if(rowOrCol == ROW)
                nextLink = rowLinks[index-rowMin];
            else
                nextLink = colLinks[index-colMin];

            linkIndex = StsEdgeLoopRadialGridLink.getLinkDirection(rowOrCol, PLUS);
        }

        public boolean hasMoreElements()
        {
            return nextLink != null;
        }

        public Object nextElement()
        {
            StsEdgeLoopRadialGridLink link;

            if(nextLink == null)
                return null;
            else
            {
                link = nextLink;
                nextLink = link.getLink(linkIndex);
                return link;
            }
        }
    }

    public boolean addInteriorLinks(StsEdgeLoop edgeLoop)
    {
        StsEdgeLoopRadialGridLink minusLink, plusLink;
        int row, col;
        StsGridSectionPoint gridPoint;
        StsEdgeLoopRadialGridLink link;

        try
        {
            for(row = rowMin; row <= rowMax; row++)
            {
                minusLink = getRowLink(row);
                while(minusLink != null)
                {
					plusLink = minusLink.getLink(ROW_PLUS);
					if(minusLink.isConnected(ROW_PLUS))
                    {
                        int colStart = StsMath.above(minusLink.getCol());
                        int colEnd = StsMath.below(plusLink.getCol());

                        for(col = colStart; col <= colEnd; col++)
                        {
                            StsPoint point = grid.getPoint(row, col);
                            gridPoint = new StsGridSectionPoint(point, row, col, null, grid, false);
                            link = StsEdgeLoopRadialGridLink.constructInteriorLink(edgeLoop, gridPoint);
                            insertOrderedLink(link);
                        }
                    }
                    minusLink = plusLink;
                }
            }
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsLinkGrid.addInteriorLinks() failed." + "\n" +
				"Grid: " + grid.getLabel(),
                e, StsException.WARNING);
            return false;
        }
    }

/*
    public final Enumeration gridLineEnum(int rowOrCol, int rowCol)
    {
        if(rowOrCol == ROW)
        {
            if(rowCol < rowMin || rowCol > rowMax)
                return new GridLineEnumerator();
            else
                return new GridLineEnumerator(ROW, rowCol);
        }
        else if(rowOrCol == COL)
        {
            if(rowCol < colMin || rowCol > colMax)
                return new GridLineEnumerator();
            else
                return new GridLineEnumerator(COL, rowCol);
        }
        else
            return new GridLineEnumerator();
    }

    private boolean insideGridRange(int rowOrCol, float rowColF)
    {
        if(rowOrCol == ROW)
            return rowColF >= rowMin && rowColF < rowMax;
        else if(rowOrCol == COL)
            return rowColF >= colMin && rowColF < colMax;
        else
            return false;
    }

    final class GridLineEnumerator implements Enumeration
    {
        int rowOrCol;
        int rowCol;
        StsLinkGridLine nextGridLine = null;
        int linkIndex;

        GridLineEnumerator()
        {
        }

        GridLineEnumerator(int rowOrCol, int rowCol)
        {
            StsEdgeLoopRadialGridLink firstLink;

            try
            {
                if(!insideGridRange(rowOrCol, rowCol)) return;

                this.rowOrCol = rowOrCol;
                this.rowCol = rowCol;

                if(rowOrCol == ROW)
                    firstLink = rowLinks[rowCol-rowMin];
                else
                    firstLink = colLinks[rowCol-colMin];

                linkIndex = StsEdgeLoopRadialGridLink.getLinkIndex(rowOrCol, PLUS);
                nextGridLine = getNextGridLine(firstLink);
            }
            catch(Exception e)
            {
                StsException.outputException("StsLinkGrid.GridLineEnumerator constructor failed.",
                    e, StsException.WARNING);
            }
        }

        private StsLinkGridLine getNextGridLine(StsEdgeLoopRadialGridLink lastLink)
        {
            while(lastLink != null)
            {
                StsEdgeLoopRadialGridLink firstLink = lastLink;
                lastLink = firstLink.getNextOutsideLink(linkIndex);
                if(lastLink != null)
                    return new StsLinkGridLine(rowOrCol, rowCol, firstLink, lastLinksystem);
                else
                    lastLink = firstLink.getLink(linkIndex);
            }
            return null;
        }

        public boolean hasMoreElements()
        {
            return nextGridLine != null;
        }

        public Object nextElement()
        {
            if(nextGridLine == null) return null;
            StsLinkGridLine gridLine = nextGridLine;
            nextGridLine = getNextGridLine(gridLine.getLastLink());
            return gridLine;
        }
    }
*/
    public StsEdgeLoopRadialGridLink getLink(int row, int col)
    {
        float colF = (float)col;
        StsEdgeLoopRadialGridLink link = rowLinks[row-rowMin];
        while(link != null)
        {
            if(link.getColF() == colF) return link;
            link = link.getLink(ROW, PLUS);
        }
        return null;
    }

    public boolean linkExists(int row, int col)
    {
        float colF = (float)col;
        StsEdgeLoopRadialGridLink link = rowLinks[row-rowMin];
        while(link != null)
        {
            if(link.getColF() == colF) return true;
            link = link.getLink(ROW, PLUS);
        }
        return false;
    }

    public StsEdgeLoopRadialGridLink getStartLink()
    {
        for(int row = rowMin; row <= rowMax; row++)
            if(hasRow(row)) return getFirstRowLink(row);

        for(int col = colMin; col <= colMax; col++)
            if(hasCol(col)) return getFirstColLink(col);

        StsException.systemError("StsLinkGrid.getStartLink() failed." +
            " There is no startLink for this linkGrid." + "\n" +
			"Grid: " + grid.getLabel());
        return null;
    }

    public void resetMinMaxLimits()
    {
        int newRowMin = rowMax;
        int newRowMax = rowMin;

        for(int row = rowMin; row <= rowMax; row++)
        {
            if(rowLinks[row-rowMin] != null)
            {
                newRowMin = row;
                break;
            }
        }

        for(int row = rowMax; row >= rowMin; row--)
        {
            if(rowLinks[row-rowMin] != null)
            {
                newRowMax = row;
                break;
            }
        }

        int nRows = rowMax-rowMin+1;
        int nNewRows = newRowMax-newRowMin+1;
        if(nRows > 0 && nNewRows > 0)
        {
            StsEdgeLoopRadialGridLink[] newRowLinks = new StsEdgeLoopRadialGridLink[nNewRows];
            System.arraycopy(rowLinks, newRowMin, newRowLinks, 0, nNewRows);
            rowMin = newRowMin;
            rowMax = newRowMax;
            rowLinks = newRowLinks;
        }
        else
        {
            rowMin = 0;
            rowMax = 1;
            rowLinks = new StsEdgeLoopRadialGridLink[2];
        }

        int newColMin = colMax;
        int newColMax = colMin;
        for(int col = colMin; col <= colMax; col++)
        {
            if(colLinks[col-colMin] != null)
            {
                newColMin = col;
                break;
            }
        }

        for(int col = colMax; col >= colMin; col--)
        {
            if(colLinks[col-colMin] != null)
            {
                newColMax = col;
                break;
            }
        }

        int nCols = colMax-colMin+1;
        int nNewCols = newColMax-newColMin+1;
        if(nCols > 0 && nNewCols > 0)
        {
            StsEdgeLoopRadialGridLink[] newColLinks = new StsEdgeLoopRadialGridLink[nNewCols];
            System.arraycopy(colLinks, newColMin, newColLinks, 0, nNewCols);
            colMin = newColMin;
            colMax = newColMax;
            colLinks = newColLinks;
        }
        else
        {
            colMin = 0;
            colMax = 1;
            colLinks = new StsEdgeLoopRadialGridLink[2];
        }
    }

    private boolean isInside(int row, int col)
    {
        try
        {
            StsEdgeLoopRadialGridLink prevLink = getFirstRowColLink(ROW, row);
            int linkIndex = StsEdgeLoopRadialGridLink.ROW_PLUS;
            StsEdgeLoopRadialGridLink nextLink = prevLink.getLink(linkIndex);

            while(nextLink != null)
            {
                if(isBetween(ROW, (float)col, prevLink, nextLink)) return true;
                prevLink = nextLink;
                nextLink = prevLink.getLink(linkIndex);
            }
            return false;
        }
        catch(Exception e)
        {
            StsException.outputException("StsLinkGrid.isInside() failed. " +
                "For point at row: " + row + " col: " + col + "\n" +
				"Grid: " + grid.getLabel(),
                e, StsException.WARNING);
            return false;
        }
    }

    public void debugDisplay(StsModel model, StsGLPanel3d glPanel3d)
    {
        StsEdgeLoopRadialGridLink link;
        int linkIndex;
        boolean debug;
        boolean displayGaps = false;

        debug = model.getBooleanProperty("Row Links ID");
        if(debug) displayRows(model, glPanel3d, true);

        debug = model.getBooleanProperty("Row Links Seq");
        if(debug) displayRows(model, glPanel3d, false);

        debug = model.getBooleanProperty("Display Rows");
        if(debug) displayRowLinks(model, glPanel3d, displayGaps);

        debug = model.getBooleanProperty("Display Columns");
        if(debug) displayColLinks(model, glPanel3d, displayGaps);

        debug = model.getBooleanProperty("Column Links ID");
        if(debug) displayCols(model, glPanel3d, true);

        debug = model.getBooleanProperty("Column Links Seq");
        if(debug) displayCols(model, glPanel3d, false);
    }

// Debug display methods

    public void displayRows(StsModel model, StsGLPanel3d glPanel3d, boolean displayID)
    {
        try
        {
            StsSpectrum spectrum = model.getSpectrum("Basic");
            float[] v = new float[3];

            GL gl = glPanel3d.getGL();
            if(gl == null) return;

            gl.glDisable(GL.GL_LIGHTING);
            glPanel3d.setViewShift(gl, StsGraphicParameters.vertexShift);

            int linkIndex = StsEdgeLoopRadialGridLink.ROW_PLUS;
            for(int row = rowMin; row <= rowMax; row++)
            {
                int n = 0;
                StsEdgeLoopRadialGridLink link  = rowLinks[row-rowMin];
                while(link != null)
                {
                    link.debugDisplay(glPanel3d, spectrum, displayID, n++, ROW);
                    link = link.getLink(linkIndex);
                }
            }
            glPanel3d.resetViewShift(gl);
            gl.glEnable(GL.GL_LIGHTING);
        }
        catch(Exception e) {} // Ignore errors, since we are debugging
    }

    public void displayCols(StsModel model, StsGLPanel3d glPanel3d, boolean displayID)
    {
        try
        {
            StsSpectrum spectrum = model.getSpectrum("Basic");
            float[] v = new float[3];

            GL gl = glPanel3d.getGL();
            if(gl == null) return;

            gl.glDisable(GL.GL_LIGHTING);
            glPanel3d.setViewShift(gl, StsGraphicParameters.vertexShift);

            int linkIndex = StsEdgeLoopRadialGridLink.COL_PLUS;
            for(int col = colMin; col <= colMax; col++)
            {
                int n = 0;
                StsEdgeLoopRadialGridLink link = colLinks[col-colMin];
                while(link != null)
                {
                    link.debugDisplay(glPanel3d, spectrum, displayID, n++, COL);
                    link = link.getLink(linkIndex);
                }
            }
            glPanel3d.resetViewShift(gl);
            gl.glEnable(GL.GL_LIGHTING);
        }
        catch(Exception e) {} // Ignore errors, since we are debugging
    }

    public void displayRowLinks(StsModel model, StsGLPanel3d glPanel3d, boolean displayGaps)
    {
        try
        {
            StsSpectrum spectrum = model.getSpectrum("Basic");

            GL gl = glPanel3d.getGL();
            gl.glDisable(GL.GL_LIGHTING);
            if(gl == null) return;

            glPanel3d.setViewShift(gl, StsGraphicParameters.vertexShift);

            int linkIndex = StsEdgeLoopRadialGridLink.ROW_PLUS;
            for(int row = rowMin; row <= rowMax; row++)
            {
                int n = 0;
                StsEdgeLoopRadialGridLink link  = rowLinks[row-rowMin];
                while(link != null)
                {
                    link.debugDisplay(glPanel3d, spectrum, true, n++, ROW);
                    StsEdgeLoopRadialGridLink nextLink = link.getLink(linkIndex);
                    if(nextLink == null) break;
                    float colF = link.getColF();
                    int colStart = StsMath.above(colF);
                    colF = nextLink.getColF();
                    int colEnd = StsMath.below(colF);

                    for(int col = colStart; col <= colEnd; col++)
                    {
                        StsPoint point = grid.getPoint(row, col);
                        StsGLDraw.drawPoint(point.v, StsColor.WHITE, glPanel3d, 3, 3, 1.0, 0.0, 0.0);
                    }
                    link = nextLink;
                }
            }
            glPanel3d.resetViewShift(gl);
            gl.glEnable(GL.GL_LIGHTING);
        }
        catch(Exception e) { } //Ignore: we're debugging
    }

    public void displayColLinks(StsModel model, StsGLPanel3d glPanel3d, boolean displayGaps)
    {
        try
        {
            StsSpectrum spectrum = model.getSpectrum("Basic");

            GL gl = glPanel3d.getGL();
            gl.glDisable(GL.GL_LIGHTING);
            if(gl == null) return;

            glPanel3d.setViewShift(gl, StsGraphicParameters.vertexShift);

            int linkIndex = StsEdgeLoopRadialGridLink.COL_PLUS;
            for(int col = colMin; col <= colMax; col++)
            {
                int n = 0;
                StsEdgeLoopRadialGridLink link  = colLinks[col-colMin];
                while(link != null)
                {
                    link.debugDisplay(glPanel3d, spectrum, true, n++, COL);
                    StsEdgeLoopRadialGridLink nextLink = link.getLink(linkIndex);
                    if(nextLink == null) break;
                    float rowF = link.getRowF();
                    int rowStart = StsMath.above(rowF);
                    rowF = nextLink.getRowF();
                    int rowEnd = StsMath.below(rowF);

                    for(int row = rowStart; row <= rowEnd; row++)
                    {
                        StsPoint point = grid.getPoint(row, col);
                        StsGLDraw.drawPoint(point.v, StsColor.WHITE, glPanel3d, 3, 3, 1.0, 0.0, 0.0);
                    }
                    link = nextLink;
                }
            }
            glPanel3d.resetViewShift(gl);
            gl.glEnable(GL.GL_LIGHTING);
        }
        catch(Exception e) { } //Ignore: we're debugging
    }
}
