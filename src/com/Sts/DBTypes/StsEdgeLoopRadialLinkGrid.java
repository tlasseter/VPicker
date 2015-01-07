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

public class StsEdgeLoopRadialLinkGrid
{
    /*
    private float rowMinF = StsParameters.largeFloat;
    private float rowMaxF = -StsParameters.largeFloat;
    private float colMinF = StsParameters.largeFloat;
    private float colMaxF = -StsParameters.largeFloat;
    protected int rowMin;
    protected int rowMax;
    protected int colMin;
    protected int colMax;
    */
    protected StsGridBoundingBox loopBoundingBox;
    protected ArrayList<StsEdgeLoop> edgeLoops;
    protected StsEdgeLoopRadialGridLink[] rowLinks;
    protected StsEdgeLoopRadialGridLink[] colLinks;
    protected StsSurfaceGridable grid;
    protected StsSurfaceGridable rowColGrid;
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

    static public final int GRID_ROW_PLUS = StsParameters.GRID_ROW_PLUS;
    static public final int GRID_COL_PLUS = StsParameters.GRID_COL_PLUS;
    static public final int GRID_ROW_MINUS = StsParameters.GRID_ROW_MINUS;
    static public final int GRID_COL_MINUS = StsParameters.GRID_COL_MINUS;

    static public final byte SECTION = StsParameters.SECTION;
    static public final byte HORIZON = StsParameters.HORIZON;
    static public final byte SECTION_VERTICAL = StsParameters.SECTION_VERTICAL;

    static final int nullInteger = StsParameters.nullInteger;
    static final float largeFloat = StsParameters.largeFloat;
    static final float roundOff = StsParameters.roundOff;

    public StsEdgeLoopRadialLinkGrid()
    {
    }

    public StsEdgeLoopRadialLinkGrid(StsModelSurface surface)
    {
        this(surface, HORIZON);
    }

    private StsEdgeLoopRadialLinkGrid(StsSurfaceGridable grid, byte gridType)
    {
        this.grid = grid;
        this.gridType = gridType;

        loopBoundingBox = grid.getGridBoundingBox();

        rowLinks = new StsEdgeLoopRadialGridLink[loopBoundingBox.getNRows()];
        colLinks = new StsEdgeLoopRadialGridLink[loopBoundingBox.getNCols()];
    }

    public StsEdgeLoopRadialLinkGrid(StsSection section, StsEdgeLoop edgeLoop)
    {
        if(edgeLoops == null) edgeLoops = new ArrayList<StsEdgeLoop>();
        edgeLoops.add(edgeLoop);

        StsList loopLinks = edgeLoop.loopLinks;
        if(loopLinks == null) return;

        this.grid = section;
        if(section.isVertical())
        {
            gridType = SECTION_VERTICAL;
            rowColGrid = null;
        }
        else
        {
            gridType = SECTION;
            rowColGrid = grid;
        }

        loopBoundingBox = computeLoopBoundingBox(edgeLoop);
        int nRows = loopBoundingBox.getNRows();
        int nCols = loopBoundingBox.getNCols();
        rowLinks = new StsEdgeLoopRadialGridLink[nRows];
        colLinks = new StsEdgeLoopRadialGridLink[nCols];

        edgeLoop.addToLinkedGrid(this);
        // resetMinMaxLimits();
        orderLinks();
    }

    public StsGridBoundingBox computeLoopBoundingBox(StsEdgeLoop loop)
    {
        return computeLoopBoundingBox(loop, true);
    }

    public StsGridBoundingBox computeGridBoundingBox(StsEdgeLoop loop)
    {
        return computeLoopBoundingBox(loop, false);
    }

    public StsGridBoundingBox computeLoopBoundingBox(StsEdgeLoop loop, boolean isLoop)
    {
        float rowF, colF;

        StsList loopLinks = loop.loopLinks;
        if(loopLinks == null) return null;
        int nLinks = loopLinks.getSize();
        float rowMinF = largeFloat;
        float rowMaxF = -largeFloat;
        float colMinF = largeFloat;
        float colMaxF = -largeFloat;
        for(int n = 0; n < nLinks; n++)
        {
            StsEdgeLoopRadialGridLink link = (StsEdgeLoopRadialGridLink) loopLinks.getElement(n);
            // StsRowCol rowCol = link.getPoint().getRowCol(rowColGrid);
            if(isLoop)
            {
                rowF = link.getRowF();
                colF = link.getColF();
            }
            else
            {
                rowF = link.getGridRowF();
                colF = link.getGridColF();
            }
            if(rowF < rowMinF)
                rowMinF = rowF;
            if(rowF > rowMaxF)
                rowMaxF = rowF;
            if(colF < colMinF)
                colMinF = colF;
            if(colF > colMaxF)
                colMaxF = colF;
        }
        int rowMin = roundUp(rowMinF);
        int rowMax = roundDown(rowMaxF);
        if(rowMin > rowMax)
        {
            rowMin = 0;
            rowMax = 0;
        }
        int colMin = roundUp(colMinF);
        int colMax = roundDown(colMaxF);
        if(colMin > colMax)
        {
            colMin = 0;
            colMax = 0;
        }
        return new StsGridBoundingBox(rowMin, rowMax, colMin, colMax);
    }

    private int roundUp(float minF)
    {
        if(minF == largeFloat || minF == -largeFloat) return 0;
        return StsMath.ceiling(minF);
    }

    private int roundDown(float maxF)
    {
        if(maxF == largeFloat || maxF == -largeFloat) return 0;
        return StsMath.floor(maxF);
    }
    /*
    public StsEdgeLoopRadialLinkGrid(StsSurfaceGridable grid, StsList edgeLinks)
    {
		this(grid);

        int nLinks = edgeLinks.getSize();
        for(int n = 0; n < nLinks; n++)
        {
            StsEdgeLoopRadialGridLink link = (StsEdgeLoopRadialGridLink)edgeLinks.getElement(n);
            addLink(link);
        }
    }
    */
// Accessors

    public StsSurfaceGridable getGrid()
    {
        return grid;
    }

    public int getRowMin()
    {
        return loopBoundingBox.rowMin;
    }

    public int getRowMax()
    {
        return loopBoundingBox.rowMax;
    }

    public int getColMin()
    {
        return loopBoundingBox.colMin;
    }

    public int getColMax()
    {
        return loopBoundingBox.colMax;
    }

    public void addEdgeLoop(StsEdgeLoop edgeLoop)
    {
        if(edgeLoops == null) edgeLoops = new ArrayList<StsEdgeLoop>();
        edgeLoops.add(edgeLoop);
        edgeLoop.addToLinkedGrid(this);
    }

    public void addLink(StsEdgeLoopRadialGridLink link)
    {
        addLink(link, true);
    }

    public void addLink(StsEdgeLoopRadialGridLink link, boolean checkConnect)
    {
        StsRowCol rowColPoint = link.getRowCol();
        int rowOrCol = rowColPoint.rowOrCol;
        int rowCol = rowColPoint.getRowCol(rowOrCol);
        // setRange(link);
        if(rowOrCol == NONE) return;
        if(rowOrCol == ROW || rowOrCol == ROWCOL)
        {
            int row = rowColPoint.getRowCol(ROW);
            insertLink(ROW, row, link, checkConnect);
        }
        if(rowOrCol == COL || rowOrCol == ROWCOL)
        {
            int col = rowColPoint.getRowCol(COL);
            insertLink(COL, col, link, checkConnect);
        }
    }

    /*
        private void setRange(StsEdgeLoopRadialGridLink link)
        {
            float rowF = link.getRowF();
            float colF = link.getColF();
            rowMinF = Math.min(rowMinF, rowF);
            rowMaxF = Math.max(rowMaxF, rowF);
            colMinF = Math.min(colMinF, colF);
            colMaxF = Math.max(colMaxF, colF);
        }
    */
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
/*
    public int getGridMin(int rowOrCol)
    {
        return rowOrCol == ROW ? colMin : rowMin;
    }

    public int getGridMax(int rowOrCol)
    {
        return rowOrCol == ROW ? colMax : rowMax;
    }
*/
    public boolean hasRow(int row)
    {
        if(rowLinks == null) return false;
        if(row < loopBoundingBox.rowMin || row > loopBoundingBox.rowMax) return false;
        return rowLinks[row - loopBoundingBox.rowMin] != null;
    }

    public boolean hasCol(int col)
    {
        if(colLinks == null) return false;
        if(col < loopBoundingBox.colMin || col > loopBoundingBox.colMax) return false;
        return colLinks[col - loopBoundingBox.colMin] != null;
    }

    public StsEdgeLoopRadialGridLink getFirstRowLink(int row)
    {
        if(rowLinks == null) return null;
        if(row < loopBoundingBox.rowMin || row > loopBoundingBox.rowMax) return null;
        return rowLinks[row - loopBoundingBox.rowMin];
    }

    public StsEdgeLoopRadialGridLink getFirstColLink(int col)
    {
        if(colLinks == null) return null;
        if(col < loopBoundingBox.colMin || col > loopBoundingBox.colMax) return null;
        return colLinks[col - loopBoundingBox.colMin];
    }

    public StsEdgeLoopRadialGridLink getFirstRowColLink(int rowOrCol, int rowCol)
    {
        try
        {
            if(rowOrCol == ROW)
                return rowLinks[rowCol - loopBoundingBox.rowMin];
            else if(rowOrCol == COL)
                return colLinks[rowCol - loopBoundingBox.colMin];
            else
                return null;
        }
        catch(Exception e)
        {
            StsException.outputException("StsEdgeLoopRadialLinkGrid.getFirstRowColLink() failed." + "\n" +
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
            rowLinks[rowCol - loopBoundingBox.rowMin] = link;
        else if(rowOrCol == COL)
            colLinks[rowCol - loopBoundingBox.colMin] = link;
    }

    public void setFirstRowColLink(int rowOrCol, StsEdgeLoopRadialGridLink link)
    {
        try
        {
            if(rowOrCol == ROW)
            {
//rowColFix                int row = StsMath.roundOffInteger(link.getRowF());
                rowLinks[link.getRow() - loopBoundingBox.rowMin] = link;
            }
            else if(rowOrCol == COL)
            {
//rowColFix                int col = StsMath.roundOffInteger(link.getColF());
                colLinks[link.getCol() - loopBoundingBox.colMin] = link;
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsEdgeLoopRadialLinkGrid.setFirstRowColLink() failed." + "\n" +
                    "Grid: " + grid.getLabel(),
                    e, StsException.WARNING);
        }
    }

    public void setStartLink(StsEdgeLoopRadialGridLink link)
    {
        int rowOrCol = link.getRowOrCol();
        if(rowOrCol == ROW)
        {
            int row = (int) link.getRowF();
            rowLinks[row - loopBoundingBox.rowMin] = link;
            return;
        }
        else if(rowOrCol == COL)
        {
            int col = (int) link.getColF();
            colLinks[col - loopBoundingBox.colMin] = link;
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
    public void orderLinks()
    {
        for(int row = loopBoundingBox.rowMin; row <= loopBoundingBox.rowMax; row++)
            if(hasRow(row)) orderLinksNew(ROW, row);

        for(int col = loopBoundingBox.colMin; col <= loopBoundingBox.colMax; col++)
            if(hasCol(col)) orderLinksNew(COL, col);
    }

    private void orderLinksNew(int rowOrCol, int rowCol)
    {
        // create loopLinkSets
        ArrayList<LoopLinkSet> loopLinkSets = new ArrayList<LoopLinkSet>();
        for(StsEdgeLoop edgeLoop : edgeLoops)
            loopLinkSets.add(new LoopLinkSet(edgeLoop, rowOrCol, rowCol));

        Iterator iterator = getIterator(rowOrCol, rowCol);
        while (iterator.hasNext())
        {
            StsEdgeLoopRadialGridLink link = (StsEdgeLoopRadialGridLink) iterator.next();
            for(LoopLinkSet loopLinkSet : loopLinkSets)
            {
                if(link.getEdgeLoop() == loopLinkSet.edgeLoop)
                {
                    loopLinkSet.addLink(link);
                    break;
                }
            }
        }
        // sort the linkSets by the min crossing col or row number
        Comparator<LoopLinkSet> linkSetComparator = new LinkSetComparator(rowOrCol);
        // for each linkSet, sort into plus, mid, and minus links ordered by crossing col or row number
        // connect all the links together
        Collections.sort(loopLinkSets, linkSetComparator);
        StsEdgeLoopRadialGridLink link = null;
        for(LoopLinkSet loopLinkSet : loopLinkSets)
            link = loopLinkSet.sortLinks(link);
    }

    public StsGridBoundingBox getLoopBoundingBox()
    {
        return loopBoundingBox;
    }

    class LoopLinkSet
    {
        /**
         * edgeLoop this linkSet belongs to
         */
        StsEdgeLoop edgeLoop;
        /**
         * flag indicating whether row (1) or col (2)
         */
        int rowOrCol;
        /**
         * row or col number
         */
        int rowCol;
        /**
         * minimum rowCol number for the positive links of this set
         */
        float minRowColF = StsParameters.largeFloat;

        ArrayList<StsEdgeLoopRadialGridLink> minusLinks = new ArrayList<StsEdgeLoopRadialGridLink>();
        ArrayList<StsEdgeLoopRadialGridLink> plusLinks = new ArrayList<StsEdgeLoopRadialGridLink>();
        ArrayList<StsEdgeLoopRadialGridLink> midLinks = new ArrayList<StsEdgeLoopRadialGridLink>();

        LoopLinkSet(StsEdgeLoop edgeLoop, int rowOrCol, int rowCol)
        {
            this.edgeLoop = edgeLoop;
            this.rowOrCol = rowOrCol;
            this.rowCol = rowCol;
        }

        void addLink(StsEdgeLoopRadialGridLink link)
        {
            int direction = link.getInsideDirection(rowOrCol);
            if(direction == MINUS)
                minusLinks.add(link);
            else if(direction == PLUS)
            {
                plusLinks.add(link);
                minRowColF = Math.min(minRowColF, link.getCrossingRowColF(rowOrCol));
            }
            else
                midLinks.add(link);
        }

        StsEdgeLoopRadialGridLink sortLinks(StsEdgeLoopRadialGridLink link)
        {
            Comparator comparator = new LinkComparatorNew(rowOrCol);
            Collections.sort(minusLinks, comparator);
            Collections.sort(plusLinks, comparator);
            Collections.sort(midLinks, comparator);

            int nMinusLinks = minusLinks.size();
            int nPlusLinks = plusLinks.size();
            int nMidLinks;

            if(nMinusLinks != nPlusLinks)
            {
                StsException.systemError("StsEdgeLoopRadialLinkGrid.orderLinks() failed." +
                        " Number of minus and plus links don't compare for: " +
                        StsParameters.rowCol(rowOrCol) + " " + rowCol + "\n" +
                        "Grid: " + grid.getLabel());
                return link;
            }

            StsEdgeLoopRadialGridLink minusLink, plusLink;

            for(int p = 0; p < nPlusLinks; p++)
            {
                plusLink = plusLinks.get(p);

                // get first minusLink; we assume this must be paired with the plusLink
                minusLink = null;
                Iterator minusLinkIter = minusLinks.iterator();
                while (minusLinkIter.hasNext())
                {
                    StsEdgeLoopRadialGridLink mLink = (StsEdgeLoopRadialGridLink) minusLinkIter.next();
                    minusLink = mLink;
                    minusLinks.remove(mLink);
                    break;
                }

                if(minusLink == null)
                {
                    StsException.systemError("StsEdgeLoopRadialLinkGrid.orderLinks() failed." +
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

                StsEdgeLoopRadialGridLink[] midLinksArray = midLinks.toArray(new StsEdgeLoopRadialGridLink[0]);
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
            return link;
        }
    }

    public final class LinkSetComparator implements Comparator<LoopLinkSet>
    {
        int crossingRowOrCol;

        LinkSetComparator(int rowOrCol)
        {
            crossingRowOrCol = StsParameters.getCrossingRowOrCol(rowOrCol);
        }

        public int compare(LoopLinkSet p1, LoopLinkSet p2)
        {
            return Float.compare(p1.minRowColF, p2.minRowColF);
        }
    }

    private void orderLinks(int rowOrCol, int rowCol)
    {
        StsEdgeLoopRadialGridLink link;

        ArrayList minusLinks = new ArrayList();
        ArrayList plusLinks = new ArrayList();
        ArrayList midLinks = new ArrayList();

        Iterator iterator = getIterator(rowOrCol, rowCol);
        while (iterator.hasNext())
        {
            link = (StsEdgeLoopRadialGridLink) iterator.next();
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
        int nMidLinks;

        if(nMinusLinks != nPlusLinks)
        {
            StsException.systemError("StsEdgeLoopRadialLinkGrid.orderLinks() failed." +
                    " Number of minus and plus links don't compare for: " +
                    StsParameters.rowCol(rowOrCol) + " " + rowCol + "\n" +
                    "Grid: " + grid.getLabel());
            return;
        }

        StsEdgeLoopRadialGridLink minusLink, plusLink;
        link = null;

        for(int p = 0; p < nPlusLinks; p++)
        {
            plusLink = (StsEdgeLoopRadialGridLink) plusLinks.get(p);

            // get first minusLink which matches this edgeLoop
            minusLink = null;
            Iterator minusLinkIter = minusLinks.iterator();
            while (minusLinkIter.hasNext())
            {
                StsEdgeLoopRadialGridLink mLink = (StsEdgeLoopRadialGridLink) minusLinkIter.next();
                if(mLink.isSameEdgeLoop(plusLink))
                {
                    minusLink = mLink;
                    minusLinks.remove(mLink);
                    break;
                }
            }

            if(minusLink == null)
            {
                StsException.systemError("StsEdgeLoopRadialLinkGrid.orderLinks() failed." +
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

            StsEdgeLoopRadialGridLink[] midLinksArray = (StsEdgeLoopRadialGridLink[]) midLinks.toArray(new StsEdgeLoopRadialGridLink[0]);
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

    /**
     * LinkComparator compares two links along a row or col sorting them by:
     * 1) loop they are on;
     * 2) increasing order of crossing col or row.
     * 3) z values.
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
            StsEdgeLoopRadialGridLink p1 = (StsEdgeLoopRadialGridLink) o1;
            StsEdgeLoopRadialGridLink p2 = (StsEdgeLoopRadialGridLink) o2;
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

    public final class LinkComparatorNew implements Comparator
    {
        int rowOrCol;

        LinkComparatorNew(int rowOrCol)
        {
            this.rowOrCol = rowOrCol;
        }

        public int compare(Object o1, Object o2)
        {
            StsEdgeLoopRadialGridLink p1 = (StsEdgeLoopRadialGridLink) o1;
            StsEdgeLoopRadialGridLink p2 = (StsEdgeLoopRadialGridLink) o2;
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
    private void insertLink(int rowOrCol, int rowCol, StsEdgeLoopRadialGridLink newLink, boolean checkConnect)
    {
        if(rowOrCol == NONE) return;
        if(checkConnect && !newLink.needsConnect(rowOrCol)) return;

        StsEdgeLoopRadialGridLink link = getFirstRowColLink(rowOrCol, rowCol);
        setFirstRowColLink(rowOrCol, newLink);
        if(link == null) return;
        newLink.initialSetLink(rowOrCol, PLUS, link);
    }

    /**
     * Add a link in the corresponding ROW and or COL.  Reorder when all links added.
     * The node is at the intersection of row and col grid lines (BOTH), or on
     * a row between col lines (ROW) or on a col between row lines (COL).
     * If BOTH, it is inserted into both row and col link lists.  If ROW or COL,
     * it is inserted only into the corresponding list.
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
        if(row >= loopBoundingBox.rowMin && row <= loopBoundingBox.rowMax)
            return rowLinks[row - loopBoundingBox.rowMin];
        else
            return null;
    }

    public StsEdgeLoopRadialGridLink getColLink(int col)
    {
        if(colLinks == null) return null;
        if(col >= loopBoundingBox.colMin && col <= loopBoundingBox.colMax)
            return colLinks[col - loopBoundingBox.colMin];
        else
            return null;
    }

    public int getColMin(int row)
    {
        if(row < loopBoundingBox.rowMin || row > loopBoundingBox.rowMax) return nullInteger;

        StsEdgeLoopRadialGridLink link = rowLinks[row - loopBoundingBox.rowMin];
        if(link == null) return nullInteger;
        return StsMath.ceiling(link.getColF());
    }

    public int getColMax(int row)
    {
        if(row < loopBoundingBox.rowMin || row > loopBoundingBox.rowMax) return nullInteger;

        StsEdgeLoopRadialGridLink link = rowLinks[row - loopBoundingBox.rowMin];
        if(link == null) return nullInteger;

        int linkDirection = StsEdgeLoopRadialGridLink.ROW_PLUS;
        StsEdgeLoopRadialGridLink nextLink;
        while ((nextLink = link.getLink(linkDirection)) != null)
            link = nextLink;

        return StsMath.floor(link.getColF());
    }

    public int getRowMin(int col)
    {
        if(col < loopBoundingBox.colMin || col > loopBoundingBox.colMax) return nullInteger;

        StsEdgeLoopRadialGridLink link = colLinks[col - loopBoundingBox.colMin];
        if(link == null) return nullInteger;
        return StsMath.ceiling(link.getRowF());
    }

    public int getRowMax(int col)
    {
        if(col < loopBoundingBox.colMin || col > loopBoundingBox.colMax) return nullInteger;

        StsEdgeLoopRadialGridLink link = colLinks[col - loopBoundingBox.colMin];
        if(link == null) return nullInteger;

        int linkIndex = StsEdgeLoopRadialGridLink.COL_PLUS;
        StsEdgeLoopRadialGridLink nextLink;
        while ((nextLink = link.getLink(linkIndex)) != null)
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

    /*
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
    */
    // use this method only when adding links on grid locations
    public void insertOrderedLink(int rowOrCol, StsEdgeLoopRadialGridLink newLink)
    {
        float newRowColF, prevRowColF, nextRowColF;
        int linkDirection;
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

            linkDirection = StsEdgeLoopRadialGridLink.getLinkDirection(rowOrCol, PLUS);

            prevLink = nextLink;
            nextLink = prevLink.getLink(linkDirection);

            while (nextLink != null)
            {
                if(isBetween(rowOrCol, newLink, prevLink, nextLink))
                {
                    newLink.connectToLinks(rowOrCol, prevLink, nextLink, grid);
//                    prevLink.setRowColConnectedLink(rowOrCol, newLink);
//                    newLink.setRowColConnectedLink(rowOrCol, nextLink);
                    return;
                }
                prevLink = nextLink;
                nextLink = prevLink.getLink(linkDirection);
            }
            StsException.systemError("StsEdgeLoopRadialLinkGrid.insertOrderedLink() failed. " +
                    "Could not insert link: " + newLink.toString() + "\n" +
                    "Grid: " + grid.getLabel());
        }
        catch(Exception e)
        {
            StsException.outputException("StsEdgeLoopRadialLinkGrid.insertOrderedLink() failed. " +
                    "Could not insert link: " + newLink.toString() + "\n" +
                    "Grid: " + grid.getLabel(), e, StsException.WARNING);
        }
    }

    /**
     * Returns an enumeration of the components of this StsList.
     *
     * @return an enumeration of the components of this list.
     */

    public final Iterator getIterator(int rowOrCol, int rowCol)
    {
        if(rowOrCol == ROW)
            return rowIterator(rowCol);
        else
            return colIterator(rowCol);
    }

    public final Iterator rowIterator(int row)
    {
        if(row < loopBoundingBox.rowMin || row > loopBoundingBox.rowMax)
            return new LineIterator();
        else
            return new LineIterator(ROW, row);
    }

    public final Iterator colIterator(int col)
    {
        if(col < loopBoundingBox.colMin || col > loopBoundingBox.colMax)
            return new LineIterator();
        else
            return new LineIterator(COL, col);
    }

    public final Iterator rowGridIterator(int row)
    {
        if(row < loopBoundingBox.rowMin || row > loopBoundingBox.rowMax)
            return new LineIterator();
        else
            return new LineIterator(ROW, GRID_ROW_PLUS, row);
    }

    public final Iterator colGridIterator(int col)
    {
        if(col < loopBoundingBox.colMin || col > loopBoundingBox.colMax)
            return new LineIterator();
        else
            return new LineIterator(COL, GRID_COL_PLUS, col);
    }

    final class LineIterator implements Iterator
    {
        int index = -1;
        StsEdgeLoopRadialGridLink nextLink = null;
        int linkIndex;

        LineIterator()
        {
        }

        LineIterator(int rowOrCol, int index)
        {
            this.index = index;

            if(rowOrCol == ROW)
                nextLink = rowLinks[index - loopBoundingBox.rowMin];
            else
                nextLink = colLinks[index - loopBoundingBox.colMin];

            linkIndex = StsEdgeLoopRadialGridLink.getLinkDirection(rowOrCol, PLUS);
        }

        LineIterator(int rowOrCol, int linkIndex, int index)
        {
            this.index = index;

            if(rowOrCol == ROW)
                nextLink = rowLinks[index - loopBoundingBox.rowMin];
            else
                nextLink = colLinks[index - loopBoundingBox.colMin];

            this.linkIndex = linkIndex;
        }

        public boolean hasNext()
        {
            return nextLink != null;
        }

        public Object next()
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

        public void remove()
        {
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
            for(row = loopBoundingBox.rowMin; row <= loopBoundingBox.rowMax; row++)
            {
                minusLink = getRowLink(row);
                while (minusLink != null)
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
            StsException.outputException("StsEdgeLoopRadialLinkGrid.addInteriorLinks() failed." + "\n" +
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
                    StsException.outputException("StsEdgeLoopRadialLinkGrid.GridLineEnumerator constructor failed.",
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
        float colF = (float) col;
        StsEdgeLoopRadialGridLink link = rowLinks[row - loopBoundingBox.rowMin];
        while (link != null)
        {
            if(link.getColF() == colF) return link;
            link = link.getLink(ROW, PLUS);
        }
        return null;
    }

    public boolean linkExists(int row, int col)
    {
        float colF = (float) col;
        StsEdgeLoopRadialGridLink link = rowLinks[row - loopBoundingBox.rowMin];
        while (link != null)
        {
            if(link.getColF() == colF) return true;
            link = link.getLink(ROW, PLUS);
        }
        return false;
    }

    public StsEdgeLoopRadialGridLink getStartLink()
    {
        for(int row = loopBoundingBox.rowMin; row <= loopBoundingBox.rowMax; row++)
            if(hasRow(row)) return getFirstRowLink(row);

        for(int col = loopBoundingBox.colMin; col <= loopBoundingBox.colMax; col++)
            if(hasCol(col)) return getFirstColLink(col);

        StsException.systemError("StsEdgeLoopRadialLinkGrid.getStartLink() failed." +
                " There is no startLink for this linkGrid." + "\n" +
                "Grid: " + grid.getLabel());
        return null;
    }

    public void resetMinMaxLimits()
    {
        int newRowMax = loopBoundingBox.rowMin;

        int newRowMin = loopBoundingBox.rowMax;
        for(int row = loopBoundingBox.rowMin; row <= loopBoundingBox.rowMax; row++)
        {
            if(rowLinks[row - loopBoundingBox.rowMin] != null)
            {
                newRowMin = row;
                break;
            }
        }
        for(int col = loopBoundingBox.colMin; col <= loopBoundingBox.colMax; col++)
        {
            if(colLinks[col - loopBoundingBox.colMin] != null)
            {
                if(colLinks[col - loopBoundingBox.colMin].getRowF() < newRowMin)
                {
                    newRowMin--;
                    break;
                }
            }
        }

        for(int row = loopBoundingBox.rowMax; row >= loopBoundingBox.rowMin; row--)
        {
            if(rowLinks[row - loopBoundingBox.rowMin] != null)
            {
                newRowMax = row;
                break;
            }
        }

        int nRows = loopBoundingBox.rowMax - loopBoundingBox.rowMin + 1;
        int nNewRows = newRowMax - newRowMin + 1;
        if(nRows > 0 && nNewRows > 0)
        {
            StsEdgeLoopRadialGridLink[] newRowLinks = new StsEdgeLoopRadialGridLink[nNewRows];
            System.arraycopy(rowLinks, newRowMin, newRowLinks, 0, nNewRows);
            loopBoundingBox.rowMin = newRowMin;
            loopBoundingBox.rowMax = newRowMax;
            rowLinks = newRowLinks;
        }
        else
        {
            loopBoundingBox.rowMin = 0;
            loopBoundingBox.rowMax = 1;
            rowLinks = new StsEdgeLoopRadialGridLink[2];
        }

        int newColMin = loopBoundingBox.colMax;
        int newColMax = loopBoundingBox.colMin;
        for(int col = loopBoundingBox.colMin; col <= loopBoundingBox.colMax; col++)
        {
            if(colLinks[col - loopBoundingBox.colMin] != null)
            {
                newColMin = col;
                break;
            }
        }

        for(int col = loopBoundingBox.colMax; col >= loopBoundingBox.colMin; col--)
        {
            if(colLinks[col - loopBoundingBox.colMin] != null)
            {
                newColMax = col;
                break;
            }
        }

        int nCols = loopBoundingBox.colMax - loopBoundingBox.colMin + 1;
        int nNewCols = newColMax - newColMin + 1;
        if(nCols > 0 && nNewCols > 0)
        {
            StsEdgeLoopRadialGridLink[] newColLinks = new StsEdgeLoopRadialGridLink[nNewCols];
            System.arraycopy(colLinks, newColMin, newColLinks, 0, nNewCols);
            loopBoundingBox.colMin = newColMin;
            loopBoundingBox.colMax = newColMax;
            colLinks = newColLinks;
        }
        else
        {
            loopBoundingBox.colMin = 0;
            loopBoundingBox.colMax = 1;
            colLinks = new StsEdgeLoopRadialGridLink[2];
        }
    }

    private boolean isInside(int row, int col)
    {
        try
        {
            StsEdgeLoopRadialGridLink prevLink = getFirstRowColLink(ROW, row);
            int linkDirection = StsEdgeLoopRadialGridLink.ROW_PLUS;
            StsEdgeLoopRadialGridLink nextLink = prevLink.getLink(linkDirection);

            while (nextLink != null)
            {
                if(isBetween(ROW, (float) col, prevLink, nextLink)) return true;
                prevLink = nextLink;
                nextLink = prevLink.getLink(linkDirection);
            }
            return false;
        }
        catch(Exception e)
        {
            StsException.outputException("StsEdgeLoopRadialLinkGrid.isInside() failed. " +
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

            int linkDirection = StsEdgeLoopRadialGridLink.ROW_PLUS;
            for(int row = loopBoundingBox.rowMin; row <= loopBoundingBox.rowMax; row++)
            {
                int n = 0;
                StsEdgeLoopRadialGridLink link = rowLinks[row - loopBoundingBox.rowMin];
                while (link != null)
                {
                    link.debugDisplay(glPanel3d, spectrum, displayID, n++, ROW);
                    link = link.getLink(linkDirection);
                }
            }
            glPanel3d.resetViewShift(gl);
            gl.glEnable(GL.GL_LIGHTING);
        }
        catch(Exception e)
        {
        } // Ignore errors, since we are debugging
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

            int linkDirection = StsEdgeLoopRadialGridLink.COL_PLUS;
            for(int col = loopBoundingBox.colMin; col <= loopBoundingBox.colMax; col++)
            {
                int n = 0;
                StsEdgeLoopRadialGridLink link = colLinks[col - loopBoundingBox.colMin];
                while (link != null)
                {
                    link.debugDisplay(glPanel3d, spectrum, displayID, n++, COL);
                    link = link.getLink(linkDirection);
                }
            }
            glPanel3d.resetViewShift(gl);
            gl.glEnable(GL.GL_LIGHTING);
        }
        catch(Exception e)
        {
        } // Ignore errors, since we are debugging
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

            int linkDirection = StsEdgeLoopRadialGridLink.ROW_PLUS;
            for(int row = loopBoundingBox.rowMin; row <= loopBoundingBox.rowMax; row++)
            {
                int n = 0;
                StsEdgeLoopRadialGridLink link = rowLinks[row - loopBoundingBox.rowMin];
                while (link != null)
                {
                    link.debugDisplay(glPanel3d, spectrum, true, n++, ROW);
                    StsEdgeLoopRadialGridLink nextLink = link.getLink(linkDirection);
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
        catch(Exception e)
        {
        } //Ignore: we're debugging
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

            int linkDirection = StsEdgeLoopRadialGridLink.COL_PLUS;
            for(int col = loopBoundingBox.colMin; col <= loopBoundingBox.colMax; col++)
            {
                int n = 0;
                StsEdgeLoopRadialGridLink link = colLinks[col - loopBoundingBox.colMin];
                while (link != null)
                {
                    link.debugDisplay(glPanel3d, spectrum, true, n++, COL);
                    StsEdgeLoopRadialGridLink nextLink = link.getLink(linkDirection);
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
        catch(Exception e)
        {
        } //Ignore: we're debugging
    }
}
