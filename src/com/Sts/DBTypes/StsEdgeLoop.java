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
import javax.media.opengl.glu.*;
import java.util.*;

public class StsEdgeLoop
{
    protected StsList directedEdges;
    protected StsSurfaceGridable grid;
    private int index;

    StsList loopLinks;
    StsEdgeLoopRadialLinkGrid linkedGrid;

    private StsCellTypeGrid cellTypeGrid;

    StsList quadStrips;
    StsList polygons;
    GridArray2D rowGridLines, colGridLines;
    boolean isPlanar = false;
    float[] normal;
    boolean subPolygonsConstructed = false;
    String displayMode = null;

    /**
     * If this Edgeloop is on a row or column aligned vertical section, the topLevel polygons are GRID polygons
     * and there are no quadStrips for this edgeLoop.
     * Otherwise this an edgeLoop on a blockGrid, a non-vertical section, or a non row or col aligned section
     * and the geometry is defined by quadStrips and trimmed polygons of type SURFACE.
     * A SURFACE quad or polygon is subdivided into GRID subpolygons.  Each GRID subpolygon is further subdivided into
     * LAYER_LEFT and LAYER_RIGHT polygons on each side of the section (if it is faulted and has two sides).
     */
    byte polygonType;

    byte constructState = NOT_CONSTRUCTED;

    private int surfDisplayListNum = 0; // display list number (> 0)
    private int gridDisplayListNum = 0; // display list number (> 0)
    private int contourDisplayListNum = 0;
    public boolean propertyChanged = true;
    public boolean colorChanged = true;

    static int nextIndex = 0;

    static final int ROW = StsParameters.ROW;
    static final int COL = StsParameters.COL;
    static final int ROWCOL = StsParameters.ROWCOL;
    static final int NONE = StsParameters.NONE;

    static final int MINUS = StsParameters.MINUS;
    static final int PLUS = StsParameters.PLUS;

    static final byte CELL_EMPTY = StsParameters.CELL_EMPTY;
    static final byte CELL_FILL = StsParameters.CELL_FULL;
    static final byte CELL_EDGE = StsParameters.CELL_EDGE;

    static final byte NOT_CONSTRUCTED = 0;
    static final byte CONSTRUCTED = 1;
    static final byte CANT_CONSTRUCT = 2;

    static final int ROW_PLUS = StsEdgeLoopRadialGridLink.ROW_PLUS;
    static final int COL_PLUS = StsEdgeLoopRadialGridLink.COL_PLUS;
    static final int ROW_MINUS = StsEdgeLoopRadialGridLink.ROW_MINUS;
    static final int COL_MINUS = StsEdgeLoopRadialGridLink.COL_MINUS;

    static public final int GRID_ROW_PLUS = StsParameters.GRID_ROW_PLUS;
    static public final int GRID_COL_PLUS = StsParameters.GRID_COL_PLUS;
    static public final int GRID_ROW_MINUS = StsParameters.GRID_ROW_MINUS;
    static public final int GRID_COL_MINUS = StsParameters.GRID_COL_MINUS;

    public StsEdgeLoop()
    {
    }

    /**
     * If this section is vertical and along a row or column, the section grid is a single cell and the polygons
     * on this edgeLoop are GRID polygons.
     */
    public StsEdgeLoop(StsBlockGrid blockGrid)
    {
        this.grid = blockGrid;
        polygonType = StsPolygon.GRID;
        initializeIndex();
    }

    /**
     * If this section is vertical and along a row or column, the section grid is a single cell and the polygons
     * on this edgeLoop are GRID polygons. Otherwise, this edgeLoop is divided by section row and column lines
     * and the top-level polygons are SURFACE polygons which are further subidvided into GRID subpolygons.
     */
    public StsEdgeLoop(StsSection section)
    {
        this.grid = section;
        if (section.isRowColAligned())
            polygonType = StsPolygon.GRID;
        else
            polygonType = StsPolygon.SECTION;
        initializeIndex();
    }

    /** An edgeLoop built for a cursorSection is transient. */
    public StsEdgeLoop(StsCursorSection grid)
    {
        this.grid = grid;
        polygonType = StsPolygon.GRID;
        initializeIndex();
    }

    private void initializeIndex() { this.index = nextIndex++; }

    public boolean initialize(StsModel model)
    {
//        return true;
        return constructEdgeLinks();
    }

    public StsList getLoopLinks()
    {
        if (loopLinks == null && !constructEdgeLinks()) return null;
        return loopLinks;
    }

    public StsSurfaceGridable getGrid(){ return grid; }

    public void setNormal(float[] normal, boolean isPlanar)
    {
        this.normal = normal;
        StsMath.normalizeVectorReturnLength(normal);
        this.isPlanar = isPlanar;
    }

    public void setIsPlanar(boolean isPlanar){ this.isPlanar = isPlanar; }

    public void addEdge(StsEdgeLinkable edge, int direction)
    {
        addDirectedEdge(new StsDirectedEdge(edge, direction));
    }

    public void addDirectedEdge(StsDirectedEdge directedEdge)
    {
        if (directedEdges == null)
            directedEdges = new StsList(4, 2);
        directedEdges.add(directedEdge);
    }

    public boolean isClosed()
    {
        StsGridSectionPoint firstPoint, lastPoint;
        StsDirectedEdge firstEdge, lastEdge;

        firstEdge = (StsDirectedEdge) directedEdges.getFirst();
        if (firstEdge == null) return false;
        firstPoint = firstEdge.getFirstPoint();

        lastEdge = (StsDirectedEdge) directedEdges.getLast();
        if (lastEdge == null) return false;
        lastPoint = lastEdge.getLastPoint();

        return lastPoint == firstPoint;
    }

    public boolean delete()
    {
        deleteDirectedEdges();
        return true;
    }

    private void deleteDirectedEdges()
    {
        int nDirectedEdges = directedEdges.getSize();
        for (int n = 0; n < nDirectedEdges; n++)
        {
            StsDirectedEdge directedEdge = (StsDirectedEdge) directedEdges.getElement(n);
            directedEdge.delete();
        }
    }

    public void deleteEdge(StsEdgeLinkable edge)
    {
        int nDirectedEdges = directedEdges.getSize();
        for (int n = 0; n < nDirectedEdges; n++)
        {
            StsDirectedEdge directedEdge = (StsDirectedEdge) directedEdges.getElement(n);
            if (directedEdge.getEdge() == edge)
            {
                directedEdges.delete(directedEdge);
                return;
            }
        }
    }

    public boolean constructEdgeLinks()
    {
        StsGridSectionPoint gridPoint;
        StsEdgeLoopRadialGridLink link, prevLink, nextLink;

        try
        {
            loopLinks = new StsList(100, 100);

            Iterator loopPointIterator = getLoopPointIterator();

            while (loopPointIterator.hasNext())
            {
                gridPoint = (StsGridSectionPoint) loopPointIterator.next();
                link = new StsEdgeLoopRadialGridLink(this, gridPoint);
                loopLinks.add(link);
            }
            loopLinks.trimToSize();

            int nLinks = loopLinks.getSize();
            if (nLinks < 2) return false;

            link = (StsEdgeLoopRadialGridLink) loopLinks.getLast();
            for (int n = 0; n < nLinks; n++)
            {
                prevLink = link;
                link = (StsEdgeLoopRadialGridLink) loopLinks.getElement(n);
                prevLink.setNextLink(link);
                link.setPrevLink(prevLink);
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.constructEdgeLinks() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    public Iterator getLoopGridPointIterator(){ return new LoopGridPointIterator(); }

    public StsCellTypeGrid getCellTypeGrid()
    {
        return cellTypeGrid;
    }

    public int getIndex()
    {
        return index;
    }

    /**
     * LoopPointEnum enumerates over all points of all loop edges
     * but does not include auxiliary edges.
     * Don't repeat the point shared by adjoining edges.
     */
    final class LoopGridPointIterator implements Iterator
    {
        Iterator loopIterator;
        StsEdgeLoopRadialGridLink link;

        LoopGridPointIterator()
        {
            loopIterator = loopLinks.iterator();
        }

        public boolean hasNext()
        {
            while (loopIterator.hasNext())
            {
                link = (StsEdgeLoopRadialGridLink) loopIterator.next();
                if (link.isRowOrCol()) return true;
            }
            return false;
        }

        public Object next(){ return link; }

        public void remove(){ }
    }

    public Iterator getLoopPointIterator(){ return new LoopPointIterator(); }

    /**
     * LoopPointEnum enumerates over all points of all loop edges
     * but does not include auxiliary edges.
     * Don't repeat the point shared by adjoining edges.
     */
    final class LoopPointIterator implements Iterator
    {
        int nDirectedEdges;
        int directedEdgeIndex = -1;
        int pointIndex = -1;
        int direction = 0;
        int lastIndex = -1;
        Object[] points;
        int pointsLength;
        Object nextPoint;

        LoopPointIterator()
        {
            setDirectedEdges();
            setPointsList(0);
        }

        private void setDirectedEdges()
        {
            if (directedEdges == null) return;
            nDirectedEdges = directedEdges.getSize();
        }

        private void setPointsList(int directedEdgeIndex)
        {
            StsDirectedEdge directedEdge = null;

            if (directedEdges == null)
            {
                nextPoint = null;
                return;
            }

            directedEdge = (StsDirectedEdge) directedEdges.getElement(directedEdgeIndex);
            if (directedEdge == null)
            {
                nextPoint = null;
                return;
            }

            StsEdgeLinkable edge = directedEdge.getEdge();
            if (edge == null)
            {
                nextPoint = null;
                return;
            }

            StsList objectList = edge.getGridEdgePointsList();
            if (objectList == null)
            {
                nextPoint = null;
                return;
            }

            points = objectList.getList();

            // don't include last point  of each edge since it is repeated on next edge
            // if not a loop, however, include last point on last edge

            pointsLength = objectList.getSize();
            direction = directedEdge.getDirection();
            boolean includeLastPoint = false;

            if (pointsLength > 0)
            {
                this.directedEdgeIndex = directedEdgeIndex;
                if (direction == PLUS)
                {
                    pointIndex = 0;
                    if (includeLastPoint)
                        lastIndex = pointsLength;
                    else
                        lastIndex = pointsLength - 1;
                }
                else
                {
                    pointIndex = pointsLength - 1;
                    if (includeLastPoint)
                        lastIndex = -1;
                    else
                        lastIndex = 0;
                }

                nextPoint = points[pointIndex];
            }
        }

        /*
                private void XsetPointsList(int directedEdgeIndex)
                {
                    StsDirectedEdge directedEdge = null;

                    if(directedEdges == null)
                    {
                        nextPoint = null;
                        return;
                    }

                    while(directedEdgeIndex < nDirectedEdges)
                    {
                        directedEdge = (StsDirectedEdge)directedEdges.getElement(directedEdgeIndex);
                        if(isEdgeInLoop(directedEdge)) break;
                        directedEdgeIndex++;
                    }

                    if(directedEdge == null)
                    {
                        nextPoint = null;
                        return;
                    }
                    StsEdgeLinkable edge = directedEdge.getEdge();
                    if(edge == null)
                    {
                        nextPoint = null;
                        return;
                    }

                    StsObjectList objectList = edge.getEdgePoints();
                    if(objectList == null)
                    {
                        nextPoint = null;
                        return;
                    }

                    points = objectList.getList();

                    // don't include last point  of each edge since it is repeated on next edge
                    // if not a loop, however, include last point on last edge

                    pointsLength = objectList.getSize();
                    direction = directedEdge.getDirection();
                    boolean includeLastPoint = !nextEdgeInLoop();

                    if(pointsLength > 0)
                    {
                        this.directedEdgeIndex = directedEdgeIndex;
                        if(direction == PLUS)
                        {
                            pointIndex = 0;
                            if(includeLastPoint)
                                lastIndex = pointsLength;
                            else
                                lastIndex = pointsLength-1;
                        }
                        else
                        {
                            pointIndex = pointsLength-1;
                            if(includeLastPoint)
                                lastIndex = -1;
                            else
                                lastIndex = 0;
                        }

                        nextPoint = points[pointIndex];
                    }
                }
        */
        public boolean isEdgeInLoop(StsDirectedEdge directedEdge)
        {
            StsEdgeLinkable edge = directedEdge.getEdge();
            if (!(edge instanceof StsSurfaceEdge)) return true;
            StsSurfaceEdge surfaceEdge = (StsSurfaceEdge) edge;
            return !surfaceEdge.isAuxiliary();
        }

        public boolean nextEdgeInLoop()
        {
            StsDirectedEdge nextDirectedEdge;

            if (directedEdgeIndex < nDirectedEdges - 1)
                nextDirectedEdge = (StsDirectedEdge) directedEdges.getElement(directedEdgeIndex + 1);
            else
                nextDirectedEdge = (StsDirectedEdge) directedEdges.getFirst();

            return isEdgeInLoop(nextDirectedEdge);
        }

        public boolean hasNext()
        {
            return nextPoint != null;
        }

        public Object next()
        {
            Object point;

            if (nextPoint == null)
                return null;
            else
            {
                point = nextPoint;
                pointIndex += direction;
                if (pointIndex != lastIndex)
                    nextPoint = points[pointIndex];
                else
                {
                    directedEdgeIndex++;
                    setPointsList(directedEdgeIndex);
                }
                return point;
            }
        }

        public void remove(){ }
    }

    public EdgeLinksIterator getEdgeLinksIterator(StsEdgeLinkable edge, int direction)
    {
        return new EdgeLinksIterator(edge, direction);
    }

    /**
     * EdgeLinksEnum enumerates over all links belonging to this edge
     * proceeding in this direction around the loop.
     * The first and last points are shared with the prev and next edges.
     */
    final class EdgeLinksIterator implements Iterator<StsEdgeLoopRadialGridLink>
    {
        int direction;
        StsEdgeLoopRadialGridLink nextLink = null;
        int nextLinkIndex;
        int startLinkIndex, endLinkIndex;

        EdgeLinksIterator(StsEdgeLinkable edge, int direction)
        {
            if (directedEdges == null) return;
            if (loopLinks == null) return;

            this.direction = direction;
            // int nLoopLinks = loopLinks.getSize();

            int nDirectedEdges = directedEdges.getSize();
            nextLinkIndex = 0;
            for (int n = 0; n < nDirectedEdges; n++)
            {
                StsDirectedEdge directedEdge = (StsDirectedEdge) directedEdges.getElement(n);
                StsEdgeLinkable directedEdgeEdge = directedEdge.getEdge();
                int nEdgePoints = directedEdgeEdge.getGridEdgePointsList().getSize();
                if (directedEdgeEdge == edge)
                {
                    if (direction == PLUS)
                    {
                        startLinkIndex = nextLinkIndex;
                        endLinkIndex = nextLinkIndex + nEdgePoints - 1;
                    }
                    else
                    {
                        endLinkIndex = nextLinkIndex;
                        startLinkIndex = nextLinkIndex + nEdgePoints - 1;
                    }
                    nextLinkIndex = startLinkIndex;
                    nextLink = (StsEdgeLoopRadialGridLink)loopLinks.getElement(nextLinkIndex);
                    return;
                }
                nextLinkIndex += nEdgePoints - 1;
            }
        }

        public boolean hasNext()
        {
            return nextLink != null;
        }

        public StsEdgeLoopRadialGridLink next()
        {
            StsEdgeLoopRadialGridLink link;

            if (nextLink == null)
                return null;
            else
            {
                link = nextLink;
                if (nextLinkIndex == endLinkIndex)
                    nextLink = null;
                else
                {
                    nextLinkIndex += direction;
                    nextLink = (StsEdgeLoopRadialGridLink)loopLinks.getElement(nextLinkIndex);
                }
                return link;
            }
        }

        public void remove() { }
    }

    public StsEdgeLoopRadialGridLink[] getEdgeEndLinks(StsEdgeLinkable edge, int direction)
    {
        int startLinkIndex, endLinkIndex;

        if (directedEdges == null) return null;
        if (loopLinks == null) return null;

        int nDirectedEdges = directedEdges.getSize();
        int nextLinkIndex = 0;
        for (int n = 0; n < nDirectedEdges; n++)
        {
            StsDirectedEdge directedEdge = (StsDirectedEdge) directedEdges.getElement(n);
            StsEdgeLinkable directedEdgeEdge = directedEdge.getEdge();
            int nEdgePoints = directedEdgeEdge.getGridEdgePointsList().getSize();
            if (directedEdgeEdge == edge)
            {
                if (direction == PLUS)
                {
                    startLinkIndex = nextLinkIndex;
                    if (n == nDirectedEdges - 1)
                        endLinkIndex = 0;
                    else
                        endLinkIndex = nextLinkIndex + nEdgePoints - 1;
                }
                else
                {
                    endLinkIndex = nextLinkIndex;
                    if (n == nDirectedEdges - 1)
                        startLinkIndex = 0;
                    else
                        startLinkIndex = nextLinkIndex + nEdgePoints - 1;

                }
                StsEdgeLoopRadialGridLink startLink = (StsEdgeLoopRadialGridLink) loopLinks.getElement(startLinkIndex);
                StsEdgeLoopRadialGridLink endLink = (StsEdgeLoopRadialGridLink) loopLinks.getElement(endLinkIndex);
                return new StsEdgeLoopRadialGridLink[]{startLink, endLink};
            }
            nextLinkIndex += nEdgePoints - 1;
        }
        return null;
    }

    // Since edges are connected in a loop, we need only add the connecting point once,
    // so skip the first point of each edge.
    // If the loop is open, add the first point of the first edge and skip the first
    // point of all other edges.
    public void addToLinkedGrid(StsEdgeLoopRadialLinkGrid linkedGrid)
    {
        this.linkedGrid = linkedGrid;

        if (loopLinks == null) return;

        int nLinks = loopLinks.getSize();
        for (int n = 0; n < nLinks; n++)
        {
            StsEdgeLoopRadialGridLink link = (StsEdgeLoopRadialGridLink) loopLinks.getElement(n);
            link.addGridConnects();
            linkedGrid.addLink(link);
        }
    }

    public StsList getDirectedEdges(){ return directedEdges; }

    public StsList getEdges()
    {
        int nEdges = directedEdges.getSize();
        StsList edges = new StsList(nEdges);
        for (int n = 0; n < nEdges; n++)
        {
            StsDirectedEdge directedEdge = (StsDirectedEdge) directedEdges.getElement(n);
            edges.add(directedEdge.getEdge());
        }

        return edges;
    }

    public StsEdgeLinkable getFirstEdge()
    {
        StsDirectedEdge directedEdge = (StsDirectedEdge) directedEdges.getFirst();
        return directedEdge.getEdge();
    }

    // If this is an OPEN_LOOP, check that it is connected to adjoining loops

    /*
        public boolean isConnected()
        {
            StsDirectedEdge directedEdge;
            int direction;

            directedEdge = (StsDirectedEdge)directedEdges.getFirst();
            direction = directedEdge.getDirection();

            if(direction == PLUS && directedEdge.getEdge().getPrevEdge() == null)
                return false;
            else if(direction == MINUS && directedEdge.getEdge().getNextEdge() == null)
                return false;

            directedEdge = (StsDirectedEdge)directedEdges.getLast();
            direction = directedEdge.getDirection();

            if(direction == PLUS && directedEdge.getEdge().getNextEdge() == null)
                return false;
            else if(direction == MINUS && directedEdge.getEdge().getPrevEdge() == null)
                return false;

            return true;
        }
    */
    private boolean constructGrid(StsGridBoundingBox surfaceBoundingBox)
    {
        try
        {
            grid.checkConstructGridNormals();
            if (!constructQuadStrips(surfaceBoundingBox)) return false;
//            polygons = constructPolygons(loopLinks);
            constructState = CONSTRUCTED;
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.constructGrid() failed.", e, StsException.WARNING);
            constructState = CANT_CONSTRUCT;
            return false;
        }
    }

    private boolean isGridVerticalRibbon()
    {
        int nRows = grid.getNRows();
        int nCols = grid.getNCols();
        return nRows <= 2 && nCols <= 2;
    }

    public StsList getQuadStrips()
    {
        if (quadStrips == null) return new StsList();
        else return quadStrips;
    }

    private boolean constructQuadStrips(StsGridBoundingBox surfaceBoundingBox)
    {
        byte currentCellType, cellType;
        QuadStrip qStrip = null;
        try
        {
            if (isGridVerticalRibbon()) return true;
            cellTypeGrid = constructCellTypeGrid(surfaceBoundingBox);
            if (cellTypeGrid == null) return false;
            quadStrips = new StsList(50, 50);

            int cellRowMin = cellTypeGrid.getRowMin();
            int cellRowMax = cellTypeGrid.getRowMax();
            int cellColMin = cellTypeGrid.getColMin();
            int cellColMax = cellTypeGrid.getColMax();

            for (int row = cellRowMin; row <= cellRowMax; row++)
            {
                currentCellType = CELL_EMPTY;
                for (int col = cellColMin; col <= cellColMax; col++)
                {
                    cellType = cellTypeGrid.getCellType(row, col);
                    if (cellType == CELL_FILL)
                    {
                        if (cellType == currentCellType)
                            qStrip.lastCol = col + 1;
                        else
                        {
                            if (qStrip != null) quadStrips.add(qStrip);
                            qStrip = new QuadStrip();
                            qStrip.rowNumber = row;
                            qStrip.firstCol = col;
                            qStrip.lastCol = col + 1;
                            qStrip.cellType = cellType;
                            currentCellType = cellType;
                        }
                    }
                    else
                    {
                        if (qStrip != null) quadStrips.add(qStrip);
                        qStrip = null;
                        currentCellType = CELL_EMPTY;
                    }
                }
                if (qStrip != null) quadStrips.add(qStrip);
                qStrip = null;
                currentCellType = CELL_EMPTY;
            }
            if (quadStrips.getSize() <= 0)
                quadStrips = null;
            else
                quadStrips.trimToSize();  // shrink to actual allocation

            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.constructQuadStrips() failed", e, StsException.WARNING);
            return false;
        }
    }

    // Access is public only because blockGrid mainDebug output may need to build it.
    public StsCellTypeGrid constructCellTypeGrid(StsGridBoundingBox surfaceBoundingBox)
    {
        if (cellTypeGrid != null) return cellTypeGrid;
        if (grid == null) return null;
        return new StsCellTypeGrid(grid, this, surfaceBoundingBox);
    }

    public StsList getPolygons()
    {
        if (polygons == null) return new StsList();
        else return polygons;
    }

    private StsList constructGridPolygonsOld(StsList links)
    {
        int n;
        StsEdgeLoopRadialGridLink firstLink = null, prevLink, link, nextLink;
        StsPolygon polygon = null;

        try
        {
            int nLinks = links.getSize();
            if (nLinks < 3) return null;

            StsList polygons = new StsList(nLinks, 10);
            boolean[] hasInsideLink = new boolean[nLinks];
            boolean polygonsConstructed = false;

            // set link flags to max number of connected edges.
            // An edge link with no insideConnections will have 2 connected edges.
            // An edge link with an insideConnection will have 4 connected edges.
            // Offset by 30 to make flag unique.

            link = (StsEdgeLoopRadialGridLink) links.getLast();
            nextLink = (StsEdgeLoopRadialGridLink) links.getFirst();
            for (n = 0; n < nLinks; n++)
            {
                prevLink = link;
                link = nextLink;
                nextLink = (StsEdgeLoopRadialGridLink) links.getElement((n + 1) % nLinks);
                hasInsideLink[n] = link.hasInsideLink(prevLink, nextLink);
                if (hasInsideLink[n])
                    link.setFlag((byte) 34);
                else
                    link.setFlag((byte) 32);
            }

            link = (StsEdgeLoopRadialGridLink) links.getFirst();
            nextLink = (StsEdgeLoopRadialGridLink) links.getSecond();
            for (n = 1; n <= nLinks; n++) // loop over all links, repeating the first one
            {
                prevLink = link;
                link = nextLink;
                nextLink = (StsEdgeLoopRadialGridLink) links.getElement((n + 1) % nLinks);

                if (!edgeSegmentOK(prevLink, link))
                {
                    polygon = null;
                    continue;
                }

                if (polygon == null)
                {
                    firstLink = prevLink;
                    polygon = new StsPolygon(prevLink, StsPolygon.GRID);
                }
                polygon.addLink(link);

                if (hasInsideLink[n % nLinks])
                {
                    if (completePolygon(polygon, firstLink, prevLink, link, false))
                    {
                        if (polygon.checkOK(null)) polygons.add(polygon);
                        polygonsConstructed = true;
                    }
                    //                else
                    //                    debugCompletePolygon(polygon);

                    polygon = null;
                }
            }

            if (!polygonsConstructed)
            {
                polygon = new StsPolygon(links, StsPolygon.GRID, grid);
                polygons.add(polygon);
            }
            if (polygons.getSize() <= 0) polygons = null;
            return polygons;
        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.constructPolygons() failed.",
                e, StsException.WARNING);
            return polygons;
        }
        finally
        {
            if (polygons != null) polygons.trimToSize();
        }
    }

    private StsList constructPolygons(StsList links, byte polygonType, boolean checkCellType)
    {
        int n;
        StsEdgeLoopRadialGridLink firstLink = null, prevLink, link, nextLink;
        StsPolygon polygon = null;

        try
        {
            int nLinks = links.getSize();
            if (nLinks < 3) return null;

            StsList polygons = new StsList(nLinks, 10);
            boolean[] hasInsideLink = new boolean[nLinks];
            boolean polygonsConstructed = false;

            // set link flags to max number of connected edges.
            // An edge link with no insideConnections will have 2 connected edges.
            // An edge link with an insideConnection will have 4 connected edges.
            // Offset by 30 to make flag unique.

            link = (StsEdgeLoopRadialGridLink) links.getLast();
            link.constructLinkedDirections();
            nextLink = (StsEdgeLoopRadialGridLink) links.getFirst();
            nextLink.constructLinkedDirections();
            for (n = 0; n < nLinks; n++)
            {
                prevLink = link;
                link = nextLink;
                nextLink = (StsEdgeLoopRadialGridLink) links.getElement((n + 1) % nLinks);
                nextLink.constructLinkedDirections();
                /*
                    int nInsideLinks = link.getNumberOfInsideLinks(prevLink, nextLink);
                    byte newFlag = 0;
                    if(nInsideLinks == 0)
                        newFlag = (byte)32;
                    else
                        newFlag = (byte)(32 + 2*nInsideLinks);
                */
                hasInsideLink[n] = link.hasInsideLink(prevLink, nextLink);
                byte flag;
                if (hasInsideLink[n])
                    flag = (byte) 34;
                else
                    flag = (byte) 32;
                /*
                  if(newFlag != flag)
                  {
                      System.out.println("different");
                  }
                */
                link.setFlag(flag);
            }

            link = (StsEdgeLoopRadialGridLink) links.getFirst();
            nextLink = (StsEdgeLoopRadialGridLink) links.getSecond();
            for (n = 1; n <= nLinks; n++) // loop over all links, repeating the first one
            {
                prevLink = link;
                link = nextLink;
                nextLink = (StsEdgeLoopRadialGridLink) links.getElement((n + 1) % nLinks);

                if (!edgeSegmentOK(prevLink, link))
                {
                    polygon = null;
                    continue;
                }

                if (polygon == null)
                {
                    firstLink = prevLink;
                    polygon = new StsPolygon(prevLink, polygonType);
                }
                polygon.addLink(link);
                // if(link.insideDirection < StsEdgeLoopRadialGridLink.NO_LINK)
                if (hasInsideLink[n % nLinks])
                {
                    if (completePolygon(polygon, firstLink, prevLink, link, false))
                    {
                        // If this cell is already defined by a quad, ignore a subsequent polygon;
                        // i.e., add a polygon only to cells which have been determined to be edge cells
                        if (polygon.checkOK(null) && (!checkCellType || cellTypeIsEdgePolygon(polygon)))
                        {
                            polygon.setLowerLeftRowCol(grid, polygonType);
                            polygons.add(polygon);
                        }
                        polygonsConstructed = true;
                    }
                    //                else
                    //                    debugCompletePolygon(polygon);

                    polygon = null;
                }
            }

            if (!polygonsConstructed)
            {
                polygon = new StsPolygon(links, polygonType, grid);
                polygons.add(polygon);
            }
            if (polygons.getSize() <= 0) polygons = null;
            return polygons;
        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.constructPolygons() failed.",
                e, StsException.WARNING);
            return polygons;
        }
        finally
        {
            if (polygons != null) polygons.trimToSize();
        }
    }

    private boolean cellTypeIsPolygon(int row, int col)
    {
        if (cellTypeGrid == null) return true;
        return cellTypeGrid.getCellType(row, col) == CELL_EDGE;
    }

    private boolean cellTypeIsEdgePolygon(StsPolygon polygon)
    {
        if (cellTypeGrid == null) return true;
        int[] rowCol = polygon.getLowerLeftRowCol(grid);
        return cellTypeGrid.getCellType(rowCol) == CELL_EDGE;
    }

    // max number edges (2 or 4) is offset by 30 to make flag unique
    private boolean edgeSegmentOK(StsEdgeLoopRadialGridLink link0, StsEdgeLoopRadialGridLink link1)
    {
        byte nEdgesLeft0 = link0.getFlag();
        if (nEdgesLeft0 == 30) return false;
        byte nEdgesLeft1 = link1.getFlag();
        if (nEdgesLeft1 == 30) return false;
        if (nEdgesLeft0 > 30) link0.setFlag(--nEdgesLeft0);
        if (nEdgesLeft1 > 30) link1.setFlag(--nEdgesLeft1);
        return true;
    }

    private boolean completePolygon(StsPolygon polygon, StsEdgeLoopRadialGridLink firstLink,
                                    StsEdgeLoopRadialGridLink prevLink, StsEdgeLoopRadialGridLink link, boolean debug)
    {
        StsEdgeLoopRadialGridLink nextLink;

        try
        {
            StsList links = polygon.getLinks();
            int nLinks = links.getSize();
            while (true)
            {
                // prevLink = link;
                // link = nextLink;
                nextLink = getNextPolygonLink(polygon, prevLink, link, debug);
                if (nextLink == null) return false;
                if (!edgeSegmentOK(link, nextLink))
                {
                    /*
                    if(debug == false) StsException.systemError("StsEdgeLoop.completePolygon() failed. " +
                                            "Edge segment from link " + link.getLabel() + " to " + nextLink.getLabel() +
                                            " already used.");
                    return false;
                    */
                }
                if (nextLink == firstLink) return true;

                polygon.addLink(nextLink);
                nLinks++;

                if (nLinks > StsPolygon.MAX_POINTS)
                {
                    if (debug == false) StsException.systemError("StsEdgeLoop.completePolygon() failed. More than " +
                        StsPolygon.MAX_POINTS + " points for polygon: " + polygon.getLabel());
                    return false;
                }
                prevLink = link;
                link = nextLink;
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.completePolygon() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    /*
        private boolean completePolygon(StsPolygon polygon, boolean mainDebug)
        {
            StsEdgeLoopRadialGridLink firstLink, nextLink;

            StsList links = polygon.getLinks();
            int nLinks = links.getSize();
            if(nLinks < 2) return false;
            firstLink = (StsEdgeLoopRadialGridLink)links.getFirst();

            while(true)
            {
                nextLink = getNextPolygonLink(polygon, mainDebug);
                if(nextLink == null) return false;
                if(nextLink == firstLink) return true;
                polygon.addLink(nextLink);
                nLinks++;

                if(nLinks > StsPolygon.MAX_POINTS)
                {
                    if(mainDebug == false) StsException.systemError("StsEdgeLoop.completePolygon() failed. More than " +
                    StsPolygon.MAX_POINTS + " points for polygon: " + polygon.getLabel());
                    return false;
                }
            }
        }
    */
    // this is first called with mainDebug == false and then with mainDebug == true
    private void debugCompletePolygon(StsPolygon polygon)
    {
        StsList links = polygon.getLinks();
        int nLinks = links.getSize();
        for (int n = nLinks - 1; n >= 2; n--)
            links.delete(n);
        StsEdgeLoopRadialGridLink firstLink = (StsEdgeLoopRadialGridLink) links.getFirst();
        StsEdgeLoopRadialGridLink prevLink = firstLink;
        StsEdgeLoopRadialGridLink link = (StsEdgeLoopRadialGridLink) links.getSecond();
        completePolygon(polygon, firstLink, prevLink, link, true);
    }


    private StsEdgeLoopRadialGridLink getNextPolygonLink(StsPolygon polygon, StsEdgeLoopRadialGridLink prevLink,
                                                         StsEdgeLoopRadialGridLink link, boolean debug)
    {
        int linkIndex = link.getCounterClockwiseLinkIndex(prevLink);
        if (linkIndex == StsEdgeLoopRadialGridLink.NO_LINK)
        {
            if (!debug)
            {
                StsException.systemError("StsEdgeLoop.getNextPolygonLink() failed." +
                    " Couldn't get next CCW linkIndex at point: " + link.getRowCol().toString() +
                    " from prevPoint: " + prevLink.getRowCol().toString());
                 linkIndex = link.getCounterClockwiseLinkIndex(prevLink);
            }
            return null;
        }
        return getNextRowColPolygonLink(link, linkIndex);
    }

    /*
        private StsEdgeLoopRadialGridLink getNextPolygonLink(StsPolygon polygon, boolean mainDebug)
        {
            StsEdgeLoopRadialGridLink link = (StsEdgeLoopRadialGridLink)(polygon.getLinks().getLast());
            if(link == null) return null;
            StsEdgeLoopRadialGridLink prevLink = (StsEdgeLoopRadialGridLink)(polygon.getLinks().getSecondToLast());
            if(prevLink == null) return null;
            int linkIndex = link.getCounterClockwiseLinkIndex(prevLink);
            if(linkIndex < 0)
            {
                if(!mainDebug)
                {
                    StsException.systemError("StsEdgeLoop.getNextPolygonLink() failed." +
                        " Couldn't get next CCW linkIndex at point: " + link.getRowCol().getLabel() +
                        " from prevPoint: " + prevLink.getRowCol().getLabel());
                }
                return null;
            }
            return getNextRowColPolygonLink(link, linkIndex);
        }
    */
    private StsEdgeLoopRadialGridLink getNextRowColPolygonLink(StsEdgeLoopRadialGridLink link, int linkIndex)
    {
        StsEdgeLoopRadialGridLink nextLink;

        nextLink = link.getLink(linkIndex);
        if (nextLink == null)
        {
            StsException.systemError("StsEdgeLoop.getNextRowColPolygonLink() failed." +
                " Linked point from: " + link.toString() + " in direction " + linkIndex +
                " doesn't exist.");
            return null;
        }

        if (isPlanar) return nextLink;

        float rowF = link.getRowColF(ROW);
        float colF = link.getRowColF(COL);
        int row, col;

        switch (linkIndex)
        {
            case 0:
                row = (int) rowF;
                col = StsMath.above(colF);
                if (nextLink.getRowColF(COL) <= col) return nextLink;
                break;
            case 1:
                col = (int) colF;
                row = StsMath.above(rowF);
                if (nextLink.getRowColF(ROW) <= row) return nextLink;
                break;
            case 2:
                row = (int) rowF;
                col = StsMath.below(colF);
                if (nextLink.getRowColF(COL) >= col) return nextLink;
                break;
            case 3:
                col = (int) colF;
                row = StsMath.below(rowF);
                if (nextLink.getRowColF(ROW) >= row) return nextLink;
                break;
            default:
                return nextLink;
        }
        nextLink = addPolygonLinkToLinkedGrid(row, col);
        if (nextLink == null)
        {
            StsException.systemError("StsEdgeLoop.getNextRowColPolygonLink() failed." +
                " Linked point from: " + link.toString() + " in direction " + linkIndex +
                " doesn't exist.");
        }
        return nextLink;
    }

    private StsEdgeLoopRadialGridLink addPolygonLinkToLinkedGrid(int row, int col)
    {
        try
        {
            StsPoint point = grid.getPoint(row, col);
            StsGridSectionPoint newPoint = new StsGridSectionPoint(point, row, col, null, grid, false);
            StsEdgeLoopRadialGridLink newLink = StsEdgeLoopRadialGridLink.constructInteriorLink(this, newPoint);
            linkedGrid.insertOrderedLink(newLink);
            return newLink;
        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.addPolygonLinkToLinkedGrid() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

    private boolean constructGridLinesOld(StsBlockGrid blockGrid)
    {
        StsEdgeLoopRadialGridLink minusLink, plusLink;
        StsGridLine gridLine;

        try
        {
            StsGridBoundingBox boundingBox = blockGrid.getGridBoundingBox();
            rowGridLines = getGridArray2D(boundingBox, ROW);
            colGridLines = getGridArray2D(boundingBox, COL);

            int nLinks = loopLinks.getSize();

            for (int n = 0; n < nLinks; n++)
            {
                minusLink = (StsEdgeLoopRadialGridLink) loopLinks.getElement(n);
                plusLink = minusLink.getLastInsideConnectedLink(ROW, PLUS);
                if (plusLink != null)
                {
                    int row = minusLink.getRow();
                    gridLine = new StsGridLine(grid, ROW, row, minusLink, plusLink);
                    rowGridLines.add(gridLine, row);
                }
                plusLink = minusLink.getLastInsideConnectedLink(COL, PLUS);
                if (plusLink != null)
                {
                    int col = minusLink.getCol();
                    gridLine = new StsGridLine(grid, COL, col, minusLink, plusLink);
                    colGridLines.add(gridLine, col);
                }
            }

            rowGridLines.trimEnds();
            colGridLines.trimEnds();

            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.constructGridLines() failed.",
                e, StsException.WARNING);
            constructState = CANT_CONSTRUCT;
            return false;
        }
    }

    private boolean constructGridLines(StsBlockGrid blockGrid)
    {
        StsGridBoundingBox boundingBox = StsBlock.modelBoundingBox;
        rowGridLines = getGridArray2D(boundingBox, ROW);
        colGridLines = getGridArray2D(boundingBox, COL);
        StsBlock block = blockGrid.block;
        if (!constructGridLines(ROW, linkedGrid.rowLinks, rowGridLines, block)) return false;
        if (!constructGridLines(COL, linkedGrid.colLinks, colGridLines, block)) return false;
        return true;
    }

    private boolean constructGridLines(int rowOrCol, StsEdgeLoopRadialGridLink[] gridLinks, GridArray2D gridLines, StsBlock block)
    {
        StsEdgeLoopRadialGridLink minusLink, plusLink;
        StsGridLine gridLine;

        try
        {
            for (StsEdgeLoopRadialGridLink gridLink : gridLinks)
            {
                minusLink = gridLink;
                while (minusLink != null)
                {
                    plusLink = minusLink.getLastConnectedLink(rowOrCol, PLUS);
                    if (plusLink == null) break;
                    if (minusLink.getEdgeLoop() == this)
                    //StsBlock linkBlock = minusLink.getPoint().getVertex().getBlock();
                    //if(linkBlock == block)
                    {
                        int rowColNumber = minusLink.getRowOrColNumber(rowOrCol);
                        gridLine = new StsGridLine(grid, rowOrCol, rowColNumber, minusLink, plusLink);
                        gridLines.add(gridLine, rowColNumber);
                    }
                    minusLink = plusLink.getLink(rowOrCol, PLUS);
                }
            }
            gridLines.trimEnds();
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.constructGridLines() failed.",
                e, StsException.WARNING);
            constructState = CANT_CONSTRUCT;
            return false;
        }
    }

    private GridArray2D getGridArray2D(StsGridBoundingBox boundingBox, int rowOrCol)
    {
        GridArray2D array = new GridArray2D(boundingBox, rowOrCol);
        if (array.arrayList == null) return null;
        else return array;
    }

    final class GridArray2D
    {
        int min, max;
        int arrayLength;
        ArrayList[] arrayList = null;

        GridArray2D(StsGridBoundingBox boundingBox, int rowOrCol)
        {
            try
            {
                if (rowOrCol == ROW)
                {
                    min = boundingBox.rowMin;
                    max = boundingBox.rowMax;
                }
                else
                {
                    min = boundingBox.colMin;
                    max = boundingBox.colMax;
                }
                arrayLength = max - min + 1;
                if (arrayLength <= 0) return;
                arrayList = new ArrayList[arrayLength];
            }
            catch(Exception e)
            {
                StsException.outputWarningException(this, "constructor", e);
            }
        }

        void add(Object object, int index)
        {
            if (index < min || index > max)
            {
                StsException.systemError("StsEdgeLoop.GridArray2D.add() failed." + " index: " + index + " is out of range " + min + " to " + max);
                return;
            }

            int i = index - min;
            if (arrayList[i] == null) arrayList[i] = new ArrayList(4);
            arrayList[i].add(object);
        }

        void trimEnds()
        {
            int iFirst = 0;
            while (iFirst < arrayLength && arrayList[iFirst] == null) iFirst++;

            if (iFirst >= arrayLength)
            {
//                StsException.systemError("StsEdgeLoop.GridArray2D.trimEnds() found null array.");
                return;
            }

            int iLast = arrayLength - 1;
            while (iLast >= iFirst && arrayList[iLast] == null) iLast--;

            // arrayList has some nulls at one or both ends: trim them off
            int newLength = iLast - iFirst + 1;
            if (newLength < arrayLength)
            {
                ArrayList[] newArrayList = new ArrayList[newLength];
                System.arraycopy(arrayList, iFirst, newArrayList, 0, newLength);
                arrayList = newArrayList;
                arrayLength = newLength;
                max = min + iLast;
                min += iFirst;
            }
        }

        void sort(Comparator comparator)
        {
            for (int n = 0; n < arrayLength; n++)
                if (arrayList[n] != null) Collections.sort(arrayList[n], comparator);
        }

        Object[] getList(int index)
        {
            if (index < min || index > max) return null;
            ArrayList list = arrayList[index - min];
            if (list == null) return null;
            return list.toArray();
        }

        Iterator iterator()
        {
            return new GridArray2DIterator(arrayList);
        }
    }

    static public final class GridArray2DIterator implements Iterator
    {
        ArrayList[] arrayList;
        int length1, length2;
        Object[] nextList;
        Object next;
        int index1, index2;

        GridArray2DIterator(ArrayList[] arrayList)
        {
            this.arrayList = arrayList;
            length1 = arrayList.length;
            index1 = 0;
            setNext();
        }

        private void setNext()
        {
            for (; index1 < length1; index1++)
            {
                if (arrayList[index1] == null) continue;
                nextList = arrayList[index1].toArray();
                if (nextList == null) continue;
                length2 = nextList.length;
                index2 = 0;
                next = nextList[0];
                return;
            }
            next = null;
        }

        public boolean hasNext(){ return next != null; }

        public Object next()
        {
            Object current = next;
            index2++;
            if (index2 < length2)
                next = nextList[index2];
            else
            {
                index1++;
                setNext();
            }
            return current;
        }

        public void remove()
        {
        }
    }

    private boolean constructGridLinesAndPolygons(StsZoneSide zoneSide)
    {
        StsSection section = zoneSide.getSection();

        StsGridBoundingBox boundingBox = zoneSide.computeGridBoundingBox();
        if (section.isVertical())
        {
            int alignment = section.getRowColAlignment();
            return constructVerticalGridLines(boundingBox, alignment);
        }
        else
        {
            // add interior links for each section row-col intersection
            // if(!linkedGrid.addInteriorLinks(this)) return false;

            // polygons =  constructPolygons(loopLinks);
            return constructNonVerticalGridLines(zoneSide, boundingBox);
        }
    }

    public boolean constructCursorSection(int rowOrCol, int rowCol, StsBlockGrid grid)
    {
        if (!constructCursorSectionGridLines(rowOrCol, rowCol, grid)) return false;
        completePolygons(polygons); // this is called AFTER the subPolygons are built above since it deletes link info
        if (!Main.isDebug) deleteTransientArrays();
        constructState = CONSTRUCTED;
        return true;
    }

    private boolean constructCursorSectionGridLines(int rowOrCol, int rowCol, StsBlockGrid grid)
    {
        int rowMin, rowMax, colMin, colMax;

        try
        {
            if (rowOrCol == ROW)
            {
                rowMin = rowCol;
                rowMax = rowCol;
                colMin = grid.getColMin();
                colMax = grid.getColMax();
            }
            else
            {
                rowMin = grid.getRowMin();
                rowMax = grid.getRowMax();
                colMin = rowCol;
                colMax = rowCol;
            }

            StsGridBoundingBox boundingBox = new StsGridBoundingBox(rowMin, rowMax, colMin, colMax, false);
            return constructVerticalGridLines(boundingBox, rowOrCol);
        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.constructCursorSectionGridLines() failed.",
                e, StsException.WARNING);
            constructState = CANT_CONSTRUCT;
            return false;
        }
    }


    private boolean constructVerticalGridLines(StsGridBoundingBox boundingBox, int alignment)
    {
        GridArray2D rowGridLinks = null, colGridLinks = null;
        int row, col;

        try
        {
            boolean hasCrossingRows = alignment != ROW;
            boolean hasCrossingCols = alignment != COL;

            // gridLines are on rows crossing this vertical section
            if (hasCrossingRows) rowGridLinks = getGridArray2D(boundingBox, ROW);

            // gridLines are on cols crossing this vertical section
            if (hasCrossingCols) colGridLinks = getGridArray2D(boundingBox, COL);

            int nLinks = loopLinks.getSize();
            for (int n = 0; n < nLinks; n++)
            {
                StsEdgeLoopRadialGridLink link = (StsEdgeLoopRadialGridLink) loopLinks.getElement(n);
                StsGridRowCol gridRowCol = link.getPoint().getGridRowCol();
                int rowOrCol = gridRowCol.getRowOrCol();
                if (hasCrossingRows && (rowOrCol == ROW || rowOrCol == ROWCOL))
                {
                    row = gridRowCol.getRow();
                    rowGridLinks.add(link, row);
                }
                else if (hasCrossingCols && (rowOrCol == COL || rowOrCol == ROWCOL))
                {
                    col = gridRowCol.getCol();
                    colGridLinks.add(link, col);
                }
            }

            if (hasCrossingRows)
            {
                rowGridLines = getGridArray2D(boundingBox, ROW);

                rowGridLinks.sort(new LinkComparator(ROW));
                int minRow = rowGridLinks.min;
                int maxRow = rowGridLinks.max;
                for (row = minRow; row <= maxRow; row++)
                {
                    Object[] gridLinks = rowGridLinks.getList(row);
                    if (gridLinks == null) continue;
                    connectVerticalLinks(ROW, row, gridLinks);
                }
                rowGridLines.trimEnds();
            }

            if (hasCrossingCols)
            {
                colGridLines = getGridArray2D(boundingBox, COL);

                colGridLinks.sort(new LinkComparator(COL));
                int minCol = colGridLinks.min;
                int maxCol = colGridLinks.max;
                for (col = minCol; col <= maxCol; col++)
                {
                    Object[] gridLinks = colGridLinks.getList(col);
                    if (gridLinks == null) continue;
                    int nGridLinks = gridLinks.length;
                    connectVerticalLinks(COL, col, gridLinks);
                }
                colGridLines.trimEnds();
            }

            polygons = constructPolygons(loopLinks, StsPolygon.GRID, true);
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.constructGridLines() failed.",
                e, StsException.WARNING);
            constructState = CANT_CONSTRUCT;
            return false;
        }
    }

    private void connectVerticalLinks(int rowOrCol, int rowCol, Object[] links)
    {
        int minusDirection, plusDirection;
        StsEdgeLoopRadialGridLink link0, link1;

        if (links == null) return;
        int nLinks = links.length;

        //if(grid instanceof StsSection || rowOrCol == COL)
        {
            plusDirection = COL_PLUS;
            minusDirection = COL_MINUS;
        }
        /*
            else
            {
                plusDirection = ROW_PLUS;
                minusDirection = ROW_MINUS;
            }
        */
        int n = 0;
        while (n < nLinks - 1)
        {
            link0 = (StsEdgeLoopRadialGridLink) links[n];
            link1 = (StsEdgeLoopRadialGridLink) links[n + 1];

            connectLinks(link0, link1, plusDirection);
            connectLinks(link1, link0, minusDirection);

//            link0.recomputeLinkedDirections();
//            link1.recomputeLinkedDirections();

            constructVerticalGridLine(rowOrCol, rowCol, link0, link1);
            n += 2;
        }
    }

    private void constructVerticalGridLine(int rowOrCol, int rowCol, StsEdgeLoopRadialGridLink link0, StsEdgeLoopRadialGridLink link1)
    {
        StsGridLinkLine gridLine = new StsGridLinkLine(grid, rowOrCol, rowCol);
        gridLine.addLink(link0);
        gridLine.addLink(link1);

        if (rowOrCol == ROW)
            rowGridLines.add(gridLine, rowCol);
        else
            colGridLines.add(gridLine, rowCol);
    }

    /**
     * This is an trimmed loop on a zoneSide. We want to construct a series of row and column grid lines which cross this edgeLoop.
     * For each polygon around the outside, and for each quad cell on the interior, we insert points where row and column grid lines
     * cross the edges of each polygon and quad.  We construct links between pairs of common row and column points which cross the
     * polygon or quad cell.  In principle, we have linked together all the grid lines.  Most of these gridLines will cross from one
     * edge of the loop to another.  There may however be internal closed loops.  We construct gridLines by following them in from the edgeLoop.
     * We flag each link as we cross it; unflagged links indicate closed internal loops which we will get in the subsequent step.
     */
    private boolean constructNonVerticalGridLines(StsZoneSide zoneSide, StsGridBoundingBox boundingBox)
    {
        StsEdgeLoopRadialGridLink link, nextLink;
        int row, col;
        try
        {
            StsSection section = zoneSide.getSection();

            // addGridCrossingLinks(zoneSide.section);
            // add grid row and col links from edgeloop as needed for gridlines along the edge(s)

            // insert grid row or col points along each interior segment
            if (!addInteriorLinks(zoneSide)) return false;

            // For each quad-cell, cross-connect grid points.
            // A quadStrip defines the lower row and the first and last column of the strip.
            // So interate from first to second-to-last column.
            StsList tempLinks = new StsList(10, 10);

            int nStrips = 0;
            if (quadStrips != null) nStrips = quadStrips.getSize();

            for (int n = 0; n < nStrips; n++)
            {
                QuadStrip q = (QuadStrip) quadStrips.getElement(n);
                row = q.rowNumber;

                for (col = q.firstCol; col < q.lastCol; col++)
                    connectQuadCellGridLines(q, row, col, tempLinks);
            }

            int nPolygons = 0;
            if (polygons != null) nPolygons = polygons.getSize();

            for (int n = 0; n < nPolygons; n++)
            {
                StsPolygon polygon = (StsPolygon) polygons.getElement(n);
                connectPolygonGridLines(polygon, tempLinks);
            }

            // Construct gridLines by following in from the edgeLoop.
            // Flag each link as we cross it; unflagged links indicate
            // closed internal loops which we will get in the subsequent step.

            rowGridLines = getGridArray2D(boundingBox, ROW);
            colGridLines = getGridArray2D(boundingBox, COL);

            int nLinks = loopLinks.getSize();
            for (int n = 0; n < nLinks; n++)
            {
                link = (StsEdgeLoopRadialGridLink) loopLinks.getElement(n);
                constructGridLine(link);
            }
            // add any edge gridLines; these are edges on blockSides which are vertical on an intersected row or column plane section
            /*
            StsObjectList edges = getEdges();
            int nEdges = edges.getSize();
            // StsBlockGrid blockGrid = zoneSide.getZoneBlock().getTopGrid();
            for(int n = 0; n < nEdges; n++)
            {
                StsEdgeLinkable edge = (StsEdgeLinkable)edges.getElement(n);
                if(!(edge instanceof StsLineZone)) continue;
                StsLineZone lineZone = (StsLineZone)edge;
                StsLine line = lineZone.line;
                int rowOrCol = line.getRowOrCol(section);
                if(rowOrCol == NONE) continue;
                int rowCol = line.getRowCol(section, blockGrid);
                if(rowOrCol == ROW)
                    rowGridLines.add(lineZone, rowCol);
                else
                    colGridLines.add(lineZone, rowCol);
            }
            */
            rowGridLines.trimEnds();
            colGridLines.trimEnds();

            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.constructGridLines(zoneSide) failed.",
                e, StsException.WARNING);
            constructState = CANT_CONSTRUCT;
            return false;
        }
    }

    private boolean addInteriorLinks(StsZoneSide zoneSide)
    {
        int row, col;
        StsEdgeLoopRadialGridLink link, nextLink;
        StsSection section = zoneSide.getSection();

        int nLoopLinks = loopLinks.getSize();
        for (int n = 0; n < nLoopLinks; n++)
        {
            link = (StsEdgeLoopRadialGridLink) loopLinks.getElement(n);
            link.addEdgeGridRowAndColLinks();
        }

        int rowMin = linkedGrid.getRowMin();
        int rowMax = linkedGrid.getRowMax();

        for (row = rowMin; row <= rowMax; row++)
        {
            Iterator rowIterator = linkedGrid.rowIterator(row);
            while (rowIterator.hasNext())
            {
                link = (StsEdgeLoopRadialGridLink) rowIterator.next();
                if (link.isConnected(ROW_PLUS))
                {
                    nextLink = link.getLink(ROW_PLUS);
                    if (!insertGridLinks(ROW, row, link, nextLink, section)) return false;
                }
            }
        }

        int colMin = linkedGrid.getColMin();
        int colMax = linkedGrid.getColMax();

        for (col = colMin; col <= colMax; col++)
        {
            Iterator colIterator = linkedGrid.colIterator(col);
            while (colIterator.hasNext())
            {
                link = (StsEdgeLoopRadialGridLink) colIterator.next();
                if (link.isConnected(COL_PLUS))
                {
                    nextLink = link.getLink(COL_PLUS);
                    if (!insertGridLinks(COL, col, link, nextLink, section)) return false;
                }
            }
        }
        return true;
    }

    private void addGridCrossingLinks(StsSection section)
    {
        StsEdgeLoopRadialGridLink link, nextLink;
        StsGridSectionPoint point, nextPoint;
        StsGridSectionPoint gridPoint;
        StsList gridPoints;
        try
        {
            int nLoopLinks = loopLinks.getSize();
            StsList newLoopLinks = new StsList(nLoopLinks, 100);
            for (int n = 0; n < nLoopLinks; n++)
            {
                link = (StsEdgeLoopRadialGridLink) loopLinks.getElement(n);
                newLoopLinks.add(link);
                nextLink = link.getNextEdgeLink();

                point = link.getPoint();
                nextPoint = nextLink.getPoint();
                gridPoints = StsGridSectionPoint.getGridCrossings(point, nextPoint, section, false, NONE);
                int nPoints = gridPoints.getSize();

                for (int i = 0; i < nPoints; i++)
                {
                    gridPoint = (StsGridSectionPoint) gridPoints.getElement(i);
                    StsEdgeLoopRadialGridLink newLink = StsEdgeLoopRadialGridLink.constructInteriorLink(this, gridPoint);
                    linkedGrid.insertOrderedLink(newLink);
                    newLink.constructLinkedDirections();
                }
                link.constructLinkedDirections();
                nextLink.constructLinkedDirections();
            }
            loopLinks = newLoopLinks;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "addGridCrossingLinks", e);
        }
    }

    private void constructGridLine(StsEdgeLoopRadialGridLink link)
    {
        StsGridRowCol gridRowCol = link.getPoint().getGridRowCol();
        int rowOrCol = gridRowCol.getRowOrCol();
        if (rowOrCol == NONE) return;
        if (rowOrCol == ROWCOL)
        {
            constructGridLine(ROW, link, gridRowCol);
            constructGridLine(COL, link, gridRowCol);
        }
        else
            constructGridLine(rowOrCol, link, gridRowCol);
    }

    private void constructGridLine(int rowOrCol, StsEdgeLoopRadialGridLink link, StsGridRowCol gridRowCol)
    {
        StsEdgeLoopRadialGridLink nextLink, prevLink;

        nextLink = link.getInsideGridLink(rowOrCol);
        if (nextLink == null) return;

        if (!link.rowColUseOK(rowOrCol)) return;

        int rowCol = gridRowCol.getRowCol(rowOrCol);
        StsGridLinkLine gridLine = new StsGridLinkLine(grid, rowOrCol, rowCol);
        gridLine.addLink(link);
        int nLinks = 0;
        while (nextLink != null)
        {
            if (nLinks++ > 1000) return;
            gridLine.addLink(nextLink);
            prevLink = link;
            link = nextLink;
            nextLink = link.getOppositeLink(prevLink);
        }
        if (gridLine.isOK(rowOrCol))
        {
            if (rowOrCol == ROW)
                rowGridLines.add(gridLine, rowCol);
            else
                colGridLines.add(gridLine, rowCol);
        }
    }

    /*
        private StsEdgeLoopRadialGridLink getNextLink(StsEdgeLoopRadialGridLink link, int rowOrCol)
        {
            StsEdgeLoopRadialGridLink[] pointLinks = link.getPointLinks();
            if(pointLinks == null) return null;
            int nPointLinks = pointLinks.length;
            for(int n = 0; n < nPointLinks; n++)
            {
                StsEdgeLoopRadialGridLink nextLink = pointLinks[n];
                if(nextLink == null) continue;
                StsGridRowCol gridRowCol = nextLink.getPoint().getGridRowCol();
                if(gridRowCol.isRowOrCol(rowOrCol)) return nextLink;
            }
            return null;
        }
    */
    private boolean insertGridLinks(int rowOrCol, int rowCol, StsEdgeLoopRadialGridLink link, StsEdgeLoopRadialGridLink nextLink,
                                    StsSection section)
    {
        StsGridSectionPoint point, nextPoint;
        StsGridSectionPoint gridPoint;
        StsList gridPoints;

        try
        {
            point = link.getPoint();
            nextPoint = nextLink.getPoint();
            gridPoints = StsGridSectionPoint.getGridCrossings(point, nextPoint, section, false, NONE);
            int nPoints = gridPoints.getSize();

            for (int n = 0; n < nPoints; n++)
            {
                gridPoint = (StsGridSectionPoint) gridPoints.getElement(n);
                gridPoint.setRowOrColIndex(section, rowOrCol, rowCol, false);
                StsEdgeLoopRadialGridLink newLink = StsEdgeLoopRadialGridLink.constructInteriorLink(this, gridPoint);
                linkedGrid.insertOrderedLink(newLink);
                newLink.constructLinkedDirections();
            }
            link.constructLinkedDirections();
            nextLink.constructLinkedDirections();
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.insertGridLinks() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    // for this section quadCell which is inside a blockZoneSide, cross connect gridLines
    private void connectQuadCellGridLines(QuadStrip q, int row, int col, StsList tempLinks)
    {
        StsEdgeLoopRadialGridLink startLink, nextLink;
        try
        {
            startLink = linkedGrid.getLink(row, col);
            if (startLink == null) return;
            nextLink = startLink.getLink(COL_PLUS);
            if (nextLink == null) return;
            StsList polygons = connectGridLines(startLink, nextLink, tempLinks);
            q.addPolygons(polygons, col);
        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.connectQuadCellGridLines() failed.",
                e, StsException.WARNING);
        }
    }

    private void connectPolygonGridLines(StsPolygon polygon, StsList tempLinks)
    {
        StsEdgeLoopRadialGridLink startLink, nextLink;
        try
        {
            StsList links = polygon.getLinks();
            startLink = (StsEdgeLoopRadialGridLink) links.getFirst();
            if (startLink == null) return;
            nextLink = (StsEdgeLoopRadialGridLink) links.getSecond();
            StsList subPolygons = connectGridLines(startLink, nextLink, tempLinks);
            polygon.setSubPolygons(subPolygons);
        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.connectQuadCellGridLines() failed.",
                e, StsException.WARNING);
        }
    }

    // cellLinks are a scratch list, allocated once, so clear it before use.
    private StsList connectGridLines(StsEdgeLoopRadialGridLink startLink, StsEdgeLoopRadialGridLink nextLink, StsList cellLinks)
    {
        StsGridSectionPoint gridPoint;
        StsEdgeLoopRadialGridLink link, prevLink;

        try
        {
            cellLinks.clear();

            link = startLink;
            cellLinks.add(link);
            int n = 0;
            while (nextLink != null && nextLink != startLink && ++n < 1000)
            {
                cellLinks.add(nextLink);
                prevLink = link;
                link = nextLink;
                nextLink = link.getNextCounterClockwiseLink(prevLink);
            }
            return connectGridLinks(cellLinks);
        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.connectGridLines() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

    private StsList connectGridLinks(StsList boundaryLinks)
    {
        int row, col;
        int gridRowMin, gridColMin, nRows, nCols;
        StsEdgeLoopRadialGridLink link;
        StsGridRowCol gridRowCol;
        try
        {
            int nBoundaryLinks = boundaryLinks.getSize();

            StsGridBoundingBox boundingBox = new StsGridBoundingBox(false);
            for (int n = 0; n < nBoundaryLinks; n++)
            {
                link = (StsEdgeLoopRadialGridLink) boundaryLinks.getElement(n);
                gridRowCol = link.getPoint().getGridRowCol();
                if (gridRowCol.isRowOrCol()) boundingBox.addPoint(gridRowCol);
            }

            gridRowMin = boundingBox.rowMin;
            gridColMin = boundingBox.colMin;

            nRows = boundingBox.getNRows();
            nCols = boundingBox.getNCols();
            if (nRows == 0 && nCols == 0)
            {
                StsPolygon polygon = new StsPolygon(boundaryLinks, StsPolygon.GRID, grid);
                StsList polygons = new StsList(1);
                polygons.add(polygon);
                completePolygons(polygons);
                return polygons;
            }

            ArrayList[] rowLists = new ArrayList[nRows];
            for (row = 0; row < nRows; row++)
                rowLists[row] = new ArrayList();

            ArrayList[] colLists = new ArrayList[nCols];
            for (col = 0; col < nCols; col++)
                colLists[col] = new ArrayList();

            for (int n = 0; n < nBoundaryLinks; n++)
            {
                link = (StsEdgeLoopRadialGridLink) boundaryLinks.getElement(n);
                gridRowCol = link.getPoint().getGridRowCol();
                int rowOrCol = gridRowCol.getRowOrCol();
                if (rowOrCol == ROW || rowOrCol == ROWCOL)
                {
                    row = gridRowCol.getRow();
                    rowLists[row - gridRowMin].add(link);
                }
                if (rowOrCol == COL || rowOrCol == ROWCOL)
                {
                    col = gridRowCol.getCol();
                    colLists[col - gridColMin].add(link);
                }
            }

            // define grid row-col intersections
            // on each rowList, add any column intersections which exist to this array
            // on each colList, if a row intersection exists, average it the existing
            // intersection

            StsEdgeLoopRadialGridLink[][] gridIntersections = new StsEdgeLoopRadialGridLink[nRows][nCols];

            Comparator comparator;

            comparator = new LinkComparator(ROW);
            for (row = 0; row < nRows; row++)
            {
                Collections.sort(rowLists[row], comparator);
                connectLinks(ROW, rowLists[row].toArray(), gridIntersections, gridRowMin, gridColMin);
            }

            comparator = new LinkComparator(COL);
            for (col = 0; col < nCols; col++)
            {
                Collections.sort(colLists[col], comparator);
                connectLinks(COL, colLists[col].toArray(), gridIntersections, gridRowMin, gridColMin);
            }

            return constructGridSubPolygons(boundingBox, boundaryLinks, gridIntersections);

        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.connectGridLinks() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

    private StsList constructGridSubPolygons(StsGridBoundingBox boundingBox, StsList boundaryLinks, StsEdgeLoopRadialGridLink[][] gridIntersections)
    {
        int row, col, i, j;
        StsPolygon polygon;

        try
        {
            StsList polygons = constructPolygons(boundaryLinks, StsPolygon.GRID, false);

            int rowMin = boundingBox.rowMin;
            int rowMax = boundingBox.rowMax;
            int colMin = boundingBox.colMin;
            int colMax = boundingBox.colMax;

            StsEdgeLoopRadialGridLink[] links = new StsEdgeLoopRadialGridLink[4];

            for (row = rowMin, i = 0; row < rowMax; row++, i++)
            {
                for (col = colMin, j = 0; col < colMax; col++, j++)
                {
                    links[0] = gridIntersections[i][j];
                    if (links[0] == null) continue;
                    links[1] = gridIntersections[i + 1][j];
                    if (links[1] == null) continue;
                    links[2] = gridIntersections[i + 1][j + 1];
                    if (links[2] == null) continue;
                    links[3] = gridIntersections[i][j + 1];
                    if (links[3] == null) continue;
                    polygon = new StsPolygon(row, col);
                    boolean forwardOrder = isForwardOrder(links);
//                    boolean forwardOrder = links[0].getLink(GRID_COL_PLUS) == links[1];
                    polygon.addLinks(links, forwardOrder);
                    polygons.add(polygon);
                }
            }
            completePolygons(polygons); // compute rowCol relative to horizon Grid & delete links
            return polygons;
        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.constructGridSubPolygons() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

    private boolean isForwardOrder(StsEdgeLoopRadialGridLink[] links)
    {
        float[][] rowAndCols = new float[3][];
        float dR0, dC0, dR1, dC1;
        if (links == null || links.length < 3) return true;
        for (int n = 0; n < 3; n++)
            rowAndCols[n] = links[n].getRowCol().getRowAndColF();

        dR0 = rowAndCols[1][0] - rowAndCols[0][0];
        dC0 = rowAndCols[1][1] - rowAndCols[0][1];
        dR1 = rowAndCols[2][0] - rowAndCols[1][0];
        dC1 = rowAndCols[2][1] - rowAndCols[1][1];

        return dR0 * dC1 - dR1 * dC0 >= 0.0f;
    }

    public boolean constructCursorSubPolygons(StsZoneBlock zoneBlock)
    {
        QuadStrip q = null;
        float[] point, normal;
        int i = -1, j = -1;

        try
        {
            if (quadStrips != null)
            {
                int nStrips = quadStrips.getSize();
                for (int n = 0; n < nStrips; n++)
                {
                    q = (QuadStrip) quadStrips.getElement(n);
                    i = q.rowNumber;

                    if (q.polygons != null)
                        constructZoneSideLayerPolygons(q, zoneBlock, StsPolygon.LAYER_RIGHT);
                    else
                    {
                        for (j = q.firstCol; j <= q.lastCol; j++)
                            constructZoneSideLayerPolygons(q, i, j, zoneBlock, StsPolygon.LAYER_RIGHT);
                    }
                }
                subPolygonsConstructed = true;
            }

            if (polygons == null) return subPolygonsConstructed;

            int nPolygons = polygons.getSize();
            for (int n = 0; n < nPolygons; n++)
            {
                StsPolygon polygon = (StsPolygon) polygons.getElement(n);
                if (polygon != null) polygon.constructZoneSideLayerPolygons(zoneBlock, StsPolygon.LAYER_RIGHT);
            }
            subPolygonsConstructed = true;
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.constructZoneSideSubPolygons() failed." +
                "row: " + i + " col: " + j, e, StsException.WARNING);
            return false;
        }
    }

    /**
     * LinkComparator compares two links along a row or col sorting them by
     * increasing order of crossing col or row.
     * If the crossing col/row values are equal, it compares edgeLoop indexes and sorts
     * them by decreasing order of index.
     * If edgeLoops are the same or don't exist, sort by z values.
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
            try
            {
                StsEdgeLoopRadialGridLink link1 = (StsEdgeLoopRadialGridLink) o1;
                StsEdgeLoopRadialGridLink link2 = (StsEdgeLoopRadialGridLink) o2;
                StsGridSectionPoint p1 = link1.getPoint();
                StsGridSectionPoint p2 = link2.getPoint();

                if (rowOrCol == NONE || rowOrCol == ROWCOL) return compareZ(p1, p2);

                StsGridRowCol rc1 = p1.getGridRowCol();
                StsGridRowCol rc2 = p2.getGridRowCol();

                float crossingRowCol1 = rc1.getCrossingRowColF(rowOrCol);
                float crossingRowCol2 = rc2.getCrossingRowColF(rowOrCol);

                if (crossingRowCol1 > crossingRowCol2) return 1;
                else if (crossingRowCol1 < crossingRowCol2) return -1;
                else
                {
                    int comp = compareZ(p1, p2);
                    if (comp != 0) return comp;

                    StsException.systemError(this, "compare", " links are identical: " + link1.toString() + " " + link2.toString());
                    return 0;
                }
            }
            catch (Exception e)
            {
                StsException.outputException("StsEdgeLoop.compare() failed.",
                    e, StsException.WARNING);
                return 0;
            }
        }
    }

    private int compareZ(StsGridSectionPoint p1, StsGridSectionPoint p2)
    {
        float z1 = p1.getZ();
        float z2 = p2.getZ();
        if (z1 > z2)
            return 1;
        else if (z1 < z2)
            return -1;
        else
            return 0;
    }

    // Links to be connected across a polygon are connected in pairs because one point
    // enters and the next point leaves; then the third point enters and the fourth point leaves.
    // If two links are on the same row or col, they are on one of the row or column orthogonal
    // edges of the polygon; don't connect them as they already should be connected.
    private void connectLinks(int gridRowOrCol, Object[] links, StsEdgeLoopRadialGridLink[][] gridIntersections,
                              int gridRowMin, int gridColMin)
    {
        StsEdgeLoopRadialGridLink link0, link1;

        if (links == null) return;
        int nLinks = links.length;

        int n = 0;
        while (n < nLinks - 1)
        {
            link0 = (StsEdgeLoopRadialGridLink) links[n];
            link1 = (StsEdgeLoopRadialGridLink) links[n + 1];

            // skip this link if it is on same row or col as next one: its already connected
            if (link0.isOnSameRowOrCol(link1))
                n += 1;
            else
            {
                // connect these two links and add any intervening grid row & col points
                connectLinks(link0, link1, gridRowOrCol, gridIntersections, gridRowMin, gridColMin);
                n += 2;
            }
        }
    }

    private void connectLinks(StsEdgeLoopRadialGridLink link0, StsEdgeLoopRadialGridLink link1, int gridRowOrCol,
                              StsEdgeLoopRadialGridLink[][] gridIntersections, int gridRowMin, int gridColMin)
    {
        StsGridSectionPoint gridPoint0, gridPoint1, gridPoint;
        StsGridRowCol gridRowCol0, gridRowCol1, gridRowCol;
//        StsPoint point0, point1, point;
        float colMinF, colMaxF, rowMinF, rowMaxF;
        int colStart, colEnd, rowStart, rowEnd;
        int inc, direction;
        int nIntersections;
        int plusDirection, minusDirection;
        StsEdgeLoopRadialGridLink prevLink, link;
        int row, col, n;
        float f;

        try
        {
            gridPoint0 = link0.getPoint();
            gridPoint1 = link1.getPoint();
            gridRowCol0 = gridPoint0.getGridRowCol();
            gridRowCol1 = gridPoint1.getGridRowCol();
//            point0 = link0.getPoint().getPoint();
//            point1 = link1.getPoint().getPoint();

            if (gridRowOrCol == ROW)
            {
                float colF0 = gridRowCol0.getColF();
                float colF1 = gridRowCol1.getColF();
                row = gridRowCol0.getRow();

                if (colF1 > colF0)
                {
                    colStart = StsMath.above(colF0);
                    colEnd = StsMath.below(colF1);
                    nIntersections = colEnd - colStart + 1;
                    inc = 1;
                    plusDirection = GRID_ROW_PLUS;
                    minusDirection = GRID_ROW_MINUS;
//                    plusDirection = ROW_PLUS;
//                    minusDirection = ROW_MINUS;
                }
                else if (colF0 > colF1)
                {
                    colStart = StsMath.below(colF0);
                    colEnd = StsMath.above(colF1);
                    nIntersections = colStart - colEnd + 1;
//                    plusDirection = ROW_MINUS;
                    plusDirection = GRID_ROW_MINUS;
//                    minusDirection = ROW_PLUS;
                    minusDirection = GRID_ROW_PLUS;
                    inc = -1;
                }
                else
                    return;

                if (nIntersections <= 0)
                {
//                    connectLinks(link0, link1);
//                    connectLinks(link1, link0);
                    connectLinks(link0, link1, plusDirection);
                    connectLinks(link1, link0, minusDirection);
                }
                else
                {
                    for (col = colStart, n = 0; n < nIntersections; n++, col += inc)
                    {
                        f = (col - colF0) / (colF1 - colF0);
                        gridPoint = StsGridSectionPoint.sectionInterpolate(gridPoint0,
                            gridPoint1, (double) f, null, (StsSection) grid, false);
//                        point = new StsPoint();
//                        point.interpolatePoints(point0, point1, f);
//                        gridPoint = new StsGridSectionPoint(point, grid, false);
                        gridRowCol = gridPoint.getGridRowCol();
                        gridRowCol.setRowOrColIndex(ROW, row, false);
                        gridRowCol.setRowOrColIndex(COL, col, false);
                        link = StsEdgeLoopRadialGridLink.constructInteriorLink(this, gridPoint);
                        gridIntersections[row - gridRowMin][col - gridColMin] = link;
                    }

                    link = gridIntersections[row - gridRowMin][colStart - gridColMin];
//                    connectLinks(link0, link);
                    connectLinks(link0, link, plusDirection);
                    connectLinks(link, link0, minusDirection);

                    for (col = colStart + 1; col <= colEnd; col++)
                    {
                        prevLink = link;
                        link = gridIntersections[row - gridRowMin][col - gridColMin];
                        connectLinks(prevLink, link, plusDirection);
                        connectLinks(link, prevLink, minusDirection);
                    }
                    connectLinks(link, link1, plusDirection);
//                    connectLinks(link1, link);
                    connectLinks(link1, link, minusDirection);
                }
//                link0.recomputeLinkedDirections();
//                link1.recomputeLinkedDirections();
            }
            else if (gridRowOrCol == COL)
            {
                float rowF0 = gridRowCol0.getRowF();
                float rowF1 = gridRowCol1.getRowF();
                col = gridRowCol0.getCol();

                if (rowF1 > rowF0)
                {
                    rowStart = StsMath.above(rowF0);
                    rowEnd = StsMath.below(rowF1);
                    nIntersections = rowEnd - rowStart + 1;
                    inc = 1;
//                    plusDirection = COL_PLUS;
                    plusDirection = GRID_COL_PLUS;
//                    minusDirection = COL_MINUS;
                    minusDirection = GRID_COL_MINUS;
                }
                else if (rowF0 > rowF1)
                {
                    rowStart = StsMath.below(rowF0);
                    rowEnd = StsMath.above(rowF1);
                    nIntersections = rowStart - rowEnd + 1;
//                    plusDirection = COL_MINUS;
                    plusDirection = GRID_COL_MINUS;
//                    minusDirection = COL_PLUS;
                    minusDirection = GRID_COL_PLUS;
                    inc = -1;
                }
                else
                    return;

                if (nIntersections <= 0)
                {
//                    connectLinks(link0, link1);
//                    connectLinks(link1, link0);
                    connectLinks(link0, link1, plusDirection);
                    connectLinks(link1, link0, minusDirection);
                }
                else
                {
                    for (row = rowStart; row <= rowEnd; row++)
                    {
                        f = (row - rowF0) / (rowF1 - rowF0);
                        StsGridSectionPoint colGridPoint = StsGridSectionPoint.sectionInterpolate(gridPoint0,
                            gridPoint1, (double) f, null, (StsSection) grid, false);

                        // point = new StsPoint();
                        // point.interpolatePoints(point0, point1, f);
                        link = getGridIntersectionLink(gridIntersections, row - gridRowMin, col - gridColMin);
                        //                       link = gridIntersections[row-gridRowMin][col-gridColMin];
                        if (link == null)
                        {
                            StsException.systemError(this, "connectLinks", "Link is not defined at row: " + row + " col: " + col);
                            continue;
                        }
                        StsGridSectionPoint rowGridPoint = link.getPoint();
                        gridPoint = StsGridSectionPoint.sectionInterpolate(colGridPoint,
                            rowGridPoint, 0.5, null, (StsSection) grid, false);
                        gridRowCol = gridPoint.getGridRowCol();
                        gridRowCol.setRowOrColIndex(ROW, row, false);
                        gridRowCol.setRowOrColIndex(COL, col, false);
                        link.resetPoint(this, gridPoint);
                    }

                    link = getGridIntersectionLink(gridIntersections, rowStart - gridRowMin, col - gridColMin);
                    // link = gridIntersections[rowStart-gridRowMin][col-gridColMin];
                    if (link != null)
                    {
//                    connectLinks(link0, link);
                        connectLinks(link0, link, plusDirection);
                        connectLinks(link, link0, minusDirection);

                        for (row = rowStart + 1; row <= rowEnd; row++)
                        {
                            prevLink = link;
                            link = getGridIntersectionLink(gridIntersections, row - gridRowMin, col - gridColMin);
                            if (link != null)
                            {
//                            link = gridIntersections[row - gridRowMin][col - gridColMin];
                                connectLinks(prevLink, link, plusDirection);
                                connectLinks(link, prevLink, minusDirection);
                            }
                        }
                        if (link != null)
                        {
                            connectLinks(link, link1, plusDirection);
                            // connectLinks(link1, link);
                            connectLinks(link1, link, minusDirection);
                        }
                    }
                }
//                link0.constructLinkedDirections();
//                link1.constructLinkedDirections();
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.defineCombineIntersections() failed.",
                e, StsException.WARNING);
        }
    }

    private StsEdgeLoopRadialGridLink getGridIntersectionLink(StsEdgeLoopRadialGridLink[][] gridIntersections, int row, int col)
    {
        try
        {
            return gridIntersections[row][col];
        }
        catch (Exception e)
        {
            return null;
        }
    }

    /*
        private void connectLinks(StsEdgeLoopRadialGridLink link0, StsEdgeLoopRadialGridLink link1)
        {
            float dCol, dRow;
            if(link0.isRow())
            {
                dRow = link1.getRowF() - link0.getRowF();
                if(dRow > 0.0f)
                    connectLinks(link0, link1, COL_PLUS);
                else if(dRow < 0.0f)
                    connectLinks(link0, link1, COL_MINUS);
                else
                {
                    StsException.systemError("StsEdgeLoop.connectLinks(link0, link1) failed." +
                        " for edgeLoop: " + getLabel() +
                        " link0: " + link0.getLabel() + " and link1: " + link1.getLabel() +
                        " are connected along the same row.");
                    return;
                }
            }
            else if(link0.isCol())
            {
                dCol = link1.getColF() - link0.getColF();
                if(dCol > 0.0f)
                    connectLinks(link0, link1, ROW_PLUS);
                else if(dCol < 0.0f)
                    connectLinks(link0, link1, ROW_MINUS);
                else
                {
                    StsException.systemError("StsEdgeLoop.connectLinks(link0, link1) failed." +
                        " for edgeLoop: " + getLabel() +
                        " link0: " + link0.getLabel() + " and link1: " + link1.getLabel() +
                        " are connected along the same col.");
                    return;
                }
            }
            else
            {
                dRow = link1.getRowF() - link0.getRowF();
                dCol = link1.getColF() - link0.getColF();

                if(Math.abs(dRow) >= Math.abs(dCol))
                {
                    if(dRow >= 0.0f && !link0.isConnected(COL_PLUS) &&
                        connectLinks(link0, link1, COL_PLUS)) return;
                    else if(!link0.isConnected(COL_MINUS) &&
                        connectLinks(link0, link1, COL_MINUS)) return;
                    else if(dCol >= 0.0f && !link0.isConnected(ROW_PLUS) &&
                        connectLinks(link0, link1, ROW_PLUS)) return;
                    else if(!link0.isConnected(ROW_MINUS) &&
                        connectLinks(link0, link1, ROW_MINUS)) return;
                }
                else
                {
                    if(dCol >= 0.0f && !link0.isConnected(ROW_PLUS) &&
                        connectLinks(link0, link1, ROW_PLUS)) return;
                    else if(!link0.isConnected(ROW_MINUS) &&
                        connectLinks(link0, link1, ROW_MINUS)) return;
                    else if(dRow >= 0.0f && !link0.isConnected(COL_PLUS) &&
                        connectLinks(link0, link1, COL_PLUS)) return;
                    else if(!link0.isConnected(COL_MINUS) &&
                        connectLinks(link0, link1, COL_MINUS)) return;
                }

                StsException.systemError("StsEdgeLoop.connectLinks(link0, link1) failed." +
                    " for edgeLoop: " + getLabel() +
                    " link0: " + link0.getLabel() + " and link1: " + link1.getLabel() +
                    " link0 is on a row and column.");
            }
        }
    */
    private boolean connectLinks(StsEdgeLoopRadialGridLink link0, StsEdgeLoopRadialGridLink link1, int rowColIndex)
    {
        try
        {
            StsEdgeLoopRadialGridLink currentLink = link0.getLink(rowColIndex);
            if (currentLink != null)
            {
                if (currentLink == link1) return true;
                rowColIndex = StsEdgeLoopRadialGridLink.getOppositeLinkIndex(rowColIndex);
                currentLink = link0.getLink(rowColIndex);
                if (currentLink != null)
                {
                    StsException.systemError("StsEdgeLoop.connectLinks(link0, link1, rowColIndex) failed." +
                        " for edgeLoop: " + getLabel() +
                        link0.toString() + " cannot connect to: " + link1.toString() +
                        ". Already connected to: " + currentLink.toString());
                    return false;
                }
            }
            link0.setLink(rowColIndex, link1);
            link0.setConnectType(rowColIndex, StsEdgeLoopRadialGridLink.INSIDE_CONNECTED);
            /*
            if(!link0.constructLinkedDirections())
            {
                link0.setLink(rowColIndex, null);
                return false;
            }
            */
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("{StsEdgeLoop.connectLinks() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    public void subdivideXYPolygons(StsSection section)
    {
//        if(linkedGrid == null) return;
//        linkedGrid.subdivideXYPolygons(section);
    }

    public void display(StsModel model, StsGLPanel3d glPanel3d, boolean displayFill, boolean displayGrid,
                        boolean displayLayers, boolean displayProperties, StsColor color, StsPropertyVolume propertyVolume, boolean isCCW)
    {
        if (constructState != CONSTRUCTED) return;
        if (displayFill)
        {
            displaySurface(model, glPanel3d, displayLayers, displayProperties, color, propertyVolume, isCCW);
        }
        if (displayGrid) displayGrid(model, glPanel3d, displayFill, color);
    }

    /*
        public void display(StsModel model, StsBlockGrid blockGrid, boolean displayFill, boolean displayGrid, boolean displayGaps, StsColor color)
        {
    //        if(!construct(blockGrid)) return;
            GL gl = model.win3d.glPanel3d.getGL();
            if(displayFill) displaySurface(model, gl, displayGaps, color);
            if(displayGrid) displayGrid(model, gl, displayFill, color);
        }
    */
    public boolean construct(StsBlockGrid blockGrid)
    {
        if (loopLinks == null && !constructEdgeLinks()) return false;
        StsGridBoundingBox gridBoundingBox = blockGrid.getGridBoundingBox();
        if (!constructGrid(gridBoundingBox)) return false;
        polygons = constructPolygons(loopLinks, StsPolygon.SURFACE, true);
        boolean ok = constructGridLines(blockGrid);
        completePolygons(polygons);  // this is called AFTER the subPolygons are built above since it deletes link info
        if (!Main.isDebug) deleteTransientArrays();
        return ok;
    }

    // we could compute using edgeLoop.loopLinks.points, but these may not have
    // been built at this point, so use the edgeLoop.edges.points which do exist.
    public StsGridBoundingBox computeGridBoundingBoxFromEdges()
    {
        StsGridBoundingBox gridBoundingBox = new StsGridBoundingBox(false);
        computeGridBoundingBoxFromEdges(gridBoundingBox);
        return gridBoundingBox;
    }

    public void computeGridBoundingBoxFromEdges(StsGridBoundingBox gridBoundingBox)
    {
        StsGridSectionPoint gridSectionPoint;

        try
        {
            StsList edges = getEdges();
            int nEdges = edges.getSize();
            for (int e = 0; e < nEdges; e++)
            {
                StsSurfaceEdge edge = (StsSurfaceEdge) edges.getElement(e);
                StsList points = edge.getGridEdgePointsList();
                int nPoints = points.getSize();
                for (int n = 0; n < nPoints; n++)
                {
                    gridSectionPoint = (StsGridSectionPoint) points.getElement(n);
                    gridBoundingBox.addPoint(gridSectionPoint.getXYZorT());
                }
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsGridBoundingBox.compute(edgeLoop) failed.", e, StsException.WARNING);
        }
    }

    public StsGridBoundingBox computeGridBoundingBoxFromLoop()
    {
        StsGridSectionPoint gridSectionPoint;

        try
        {
            StsGridBoundingBox gridBoundingBox = new StsGridBoundingBox(false);
            int nLinks = loopLinks.getSize();
            for (int n = 0; n < nLinks; n++)
            {
                StsEdgeLoopRadialGridLink link = (StsEdgeLoopRadialGridLink) loopLinks.getElement(n);
                gridSectionPoint = link.getPoint();
                gridBoundingBox.addPoint(gridSectionPoint.getXYZorT());
            }
            return gridBoundingBox;
        }
        catch (Exception e)
        {
            StsException.outputException("StsGridBoundingBox.compute(edgeLoop) failed.", e, StsException.WARNING);
            return null;
        }
    }

    public boolean addInsidePointsToGridBoundingBox(StsGridBoundingBox gridBoundingBox, StsSurfaceGridable paramGrid, StsXYGridable xyGrid)
    {
        StsEdgeLoopRadialGridLink minusLink, plusLink;
        int row, col;

        try
        {
            int rowMin = linkedGrid.getRowMin();
            int rowMax = linkedGrid.getRowMax();
            for (row = rowMin; row <= rowMax; row++)
            {
                plusLink = linkedGrid.getRowLink(row);
                if (plusLink == null) continue;
                while (true)
                {
                    minusLink = plusLink;
                    plusLink = minusLink.getLink(ROW, PLUS);
                    if (plusLink == null) break;
                    if (minusLink.isConnected(ROW, PLUS))
                    {
                        float colMinF = minusLink.getColF();
                        int colMin = StsMath.above(colMinF);
                        float colMaxF = plusLink.getColF();
                        int colMax = StsMath.below(colMaxF);
                        for (col = colMin; col <= colMax; col++)
                        {
                            float[] xy = paramGrid.getXYZorT(row, col);
                            float gridRowF = xyGrid.getRowCoor(xy);
                            float gridColF = xyGrid.getColCoor(xy);
                            gridBoundingBox.addPoint(paramGrid.getXYZorT(gridRowF, gridColF));
                        }

                    }
                }
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsGridBoundingBox.compute(edgeLoop) failed.", e, StsException.WARNING);
            return false;
        }
    }

    public StsGridBoundingBox computeGridBoundingBox(StsSurfaceGridable paramGrid, StsXYGridable xyGrid)
    {
        StsEdgeLoopRadialGridLink minusLink, plusLink;
        StsGridSectionPoint gridSectionPoint;
        int row, col;

        try
        {
            StsGridBoundingBox gridBoundingBox = new StsGridBoundingBox(false);
            int nLinks = loopLinks.getSize();

            for (int n = 0; n < nLinks; n++)
            {
                minusLink = (StsEdgeLoopRadialGridLink) loopLinks.getElement(n);
                gridSectionPoint = minusLink.getPoint();
                gridBoundingBox.addPoint(gridSectionPoint.getXYZorT());

                plusLink = minusLink.getLastConnectedLink(ROW, PLUS);
                if (plusLink != null)
                {
                    gridSectionPoint = plusLink.getPoint();
                    gridBoundingBox.addPoint(gridSectionPoint.getXYZorT());

                    int colMin = StsMath.above(minusLink.getCol());
                    int colMax = StsMath.below(plusLink.getCol());
                    row = minusLink.getRow();

                    for (col = colMin; col <= colMax; col++)
                    {
                        gridBoundingBox.addPoint(paramGrid.getXYZorT(row, col));
                    }
                }
                plusLink = minusLink.getLastConnectedLink(COL, PLUS);
                if (plusLink != null)
                {
                    gridSectionPoint = plusLink.getPoint();
                    gridBoundingBox.addPoint(gridSectionPoint.getXYZorT());

                    int rowMin = StsMath.above(minusLink.getRow());
                    int rowMax = StsMath.below(plusLink.getRow());
                    col = minusLink.getCol();

                    for (row = rowMin; row <= rowMax; row++)
                    {
                        gridBoundingBox.addPoint(paramGrid.getXYZorT(row, col));
                    }
                }
            }
            gridBoundingBox.computeRowColRanges(xyGrid);
            return gridBoundingBox;
        }
        catch (Exception e)
        {
            StsException.outputException("StsGridBoundingBox.compute(edgeLoop) failed.", e, StsException.WARNING);
            return null;
        }
    }
/*
    public boolean construct(StsBlockGrid blockGrid)
    {
        if(!constructGrid()) return false;
        polygons = constructPolygons(loopLinks);
        boolean ok = constructGridLines(blockGrid);
        completePolygons(polygons, grid);  // this is called AFTER the subPolygons are built above since it deletes link info
        if(!Main.isDebug) deleteTransientArrays();
        return ok;
    }
*/

    /*
        public void display(StsModel model, StsZoneSide zoneSide,
                        boolean displayFill, boolean displayGrid, boolean displayGaps, StsColor color)
        {
            if(!construct(zoneSide)) return;
            GL gl = model.win3d.glPanel3d.getGL();
            if(displayFill) displaySurface(model, gl, displayGaps, color);
            if(displayGrid) displayGrid(model, gl, displayFill, color);
        }
    */
    public boolean constructPolygonsAndGridLines(StsZoneSide zoneSide)
    {

        StsGridBoundingBox surfaceBoundingBox = zoneSide.getSectionBoundingBox();
        if (!constructGrid(surfaceBoundingBox)) return false;
        if (!constructSectionPolygons(zoneSide)) return false;
        if (!constructGridLinesAndPolygons(zoneSide)) return false;
        completePolygons(polygons); // this is called AFTER the subPolygons are built above since it deletes link info
        if (!Main.isDebug) deleteTransientArrays();
        return true;
    }

    /**
     * For this zoneSide, if it is not vertical, construct polygons around trimmed edge.  Because these polygons contain
     * vertical grid intersections, they will be grid polygons even if they are not bounded by row or column grid lines.
     * For each of these polygons, assign it as type GRID and assign the row and column number as being the lower-left row & column
     * containing this polygon.
     *
     * @param zoneSide
     * @return
     */
    private boolean constructSectionPolygons(StsZoneSide zoneSide)
    {
        if (zoneSide.getSection().isVertical()) return true;
        // add interior links for each section row-col intersection
        if (!linkedGrid.addInteriorLinks(this)) return false;
        polygons = constructPolygons(loopLinks, StsPolygon.SECTION, true);
        return true;
    }

    //replace links array with xyz array in each polygon

    private void completePolygons(StsList polygons)
    {
        if (polygons == null) return;
        int nPolygons = polygons.getSize();
        for (int n = 0; n < nPolygons; n++)
        {
            StsPolygon polygon = (StsPolygon) polygons.getElement(n);
            polygon.computePoints(grid, isPlanar);
        }
    }

    /**
     * Construct subpolygons for each quad and polygon inside this trimming edgeloop on this zoneSide.
     * These subpolygons are generated by vertical grid row and column intersections of the parent quads and polygons.
     * These subpolygons are of type GRID and have a row and column number corresponding to the grid row and columns intersecting.
     *
     * @param zoneSide zoneSide defining this trimloop
     * @return
     */
    public boolean constructZoneSideLayerPolygons(StsZoneSide zoneSide)
    {
        QuadStrip q;
        int sectionRow = -1, sectionCol = -1;

        try
        {
            byte layerSide;
            if (zoneSide.side == StsParameters.RIGHT)
                layerSide = StsPolygon.LAYER_RIGHT;
            else
                layerSide = StsPolygon.LAYER_LEFT;

            StsZoneBlock zoneBlock = zoneSide.getZoneBlock();
            if (quadStrips != null)
            {
                int nStrips = quadStrips.getSize();
                for (int n = 0; n < nStrips; n++)
                {
                    q = (QuadStrip) quadStrips.getElement(n);
                    sectionRow = q.rowNumber;

                    if (q.polygons != null)
                        constructZoneSideLayerPolygons(q, zoneBlock, layerSide);
                    else
                    {
                        for (sectionCol = q.firstCol; sectionCol < q.lastCol; sectionCol++)
                            constructZoneSideLayerPolygons(q, sectionRow, sectionCol, zoneBlock, layerSide);
                    }
                }
                subPolygonsConstructed = true;
            }

            if (polygons == null) return subPolygonsConstructed;

            int nPolygons = polygons.getSize();
            for (int n = 0; n < nPolygons; n++)
            {
                StsPolygon polygon = (StsPolygon) polygons.getElement(n);
                if (polygon != null) polygon.constructZoneSideLayerPolygons(zoneBlock, layerSide);
            }
            subPolygonsConstructed = true;
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.constructZoneSideSubPolygons() failed." +
                "row: " + sectionRow + " col: " + sectionCol, e, StsException.WARNING);
            return false;
        }
    }

    private void constructZoneSideLayerPolygons(QuadStrip q, int i, int j, StsZoneBlock zoneBlock, byte layerSide)
    {
        try
        {
            float[][] pntNrmls = new float[4][];
            pntNrmls[0] = getPntNrml(i, j);
            pntNrmls[1] = getPntNrml(i + 1, j);
            pntNrmls[2] = getPntNrml(i + 1, j + 1);
            pntNrmls[3] = getPntNrml(i, j + 1);

            StsPolygon polygon = new StsPolygon();
            polygon.setPntNrmls(pntNrmls);
            polygon.constructZoneSideLayerPolygons(zoneBlock, layerSide);
            StsList polygonList = new StsList(1);
            polygonList.add(polygon);
            q.addPolygons(polygonList, j);
        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.constructZoneSideSubPolygons(q, i, j, zb) failed.",
                e, StsException.WARNING);
        }
    }

    private float[] getPntNrml(int row, int col)
    {
        float[] point = grid.getXYZorT(row, col);
        float[] normal = grid.getNormal(row, col);
        if (normal != null)
            return new float[]{point[0], point[1], point[2], normal[0], normal[1], normal[2]};
        else
            return new float[]{point[0], point[1], point[2], 0.0f, 0.0f, 1.0f};
    }

    private void constructZoneSideLayerPolygons(QuadStrip q, StsZoneBlock zoneBlock, byte layerSide)
    {
        try
        {
            if (q.polygons == null) return;

            int firstCol = q.firstCol;
            int lastCol = q.lastCol;
            int nCols = lastCol - firstCol;
            for (int j = 0; j < nCols; j++)
            {
                Object[] cellPolygons = q.polygons[j];
                if (cellPolygons == null) continue;
                int nPolygons = cellPolygons.length;
                for (int n = 0; n < nPolygons; n++)
                {
                    StsPolygon polygon = (StsPolygon) cellPolygons[n];
                    if (polygon != null) polygon.constructZoneSideLayerPolygons(zoneBlock, layerSide);
                }
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.constructZoneSideSubPolygons() failed.",
                e, StsException.WARNING);
        }
    }

    /*
        public void buildNewDomainGeometry(byte zDomain, StsSeismicVelocityModel velocityModel)
        {
            if(grid instanceof StsXYSurfaceGridable)
            {
                buildNewDomainXYSurfaceGeometry(zDomain, velocityModel, (StsXYSurfaceGridable)grid);
            }
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
                float[][] pointsZ = xyGrid.getAdjPointsZ();
                for(row = rowMin; row <= rowMax; row++, i++)
                {
                    float y = yMin + row*yInc;
                    j = 0;
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
                    "row: " + row + " col: " + col, e, StsException.WARNING);
            }
        }

        private void buildNewDomainSurfaceGeometry(byte zDomain, StsSeismicVelocityModel velocityModel)
        {
            QuadStrip q = null;
            float[] point, normal;
            int i = -1, j = -1;

            try
            {
            }
            catch(Exception e)
            {
                StsException.outputException("StsEdgeLoop.drawQuadStrips() failed." +
                    "row: " + i + " col: " + j, e, StsException.WARNING);
            }
        }
    */
    private void displaySurface(StsModel model, StsGLPanel3d glPanel3d, boolean displayLayers, boolean displayProperties, StsColor color, StsPropertyVolume propertyVolume, boolean isCCW)
    {
        GL gl = glPanel3d.getGL();
        if (gl == null) return;
        try
        {
            gl.glEnable(GL.GL_LIGHTING);  // shouldn't need this

            if (!isCCW) gl.glFrontFace(GL.GL_CW);  // change to left-handed coordinates (z is down)
            gl.glShadeModel(GL.GL_SMOOTH);
            if (StsGLPanel.debugProjectionMatrix)
                StsException.systemDebug("no viewShift called for displaySurface of " + toString());

            if(color != null) color.setGLColor(gl);

            if(propertyChanged || colorChanged)
                deleteSurfDisplayList(gl);

            boolean useDisplayLists = model.useDisplayLists;

            if (surfDisplayListNum == 0 && useDisplayLists)  // build display list
            {
                surfDisplayListNum = gl.glGenLists(1);
                if (surfDisplayListNum == 0)
                {
                    StsMessageFiles.logMessage("System Error in StsGrid.displaySurface: " +
                        "Failed to allocate a display list");
                    return;
                }

                gl.glNewList(surfDisplayListNum, GL.GL_COMPILE_AND_EXECUTE);
                drawSurface(glPanel3d, displayLayers, displayProperties, propertyVolume, propertyChanged);
                //            normals = null; // delete normals after building for display list
                gl.glEndList();
                //timer.stop("display list surface setup: ");

            }
            else if (useDisplayLists) // use existing display list
            {
                // gl.glEnable(GL.GL_LIGHTING);  // shouldn't need this
                gl.glCallList(surfDisplayListNum);
                //timer.stop("display list surface draw: ");
            }
            else // immediate mode draw
            {
                gl.glEnable(GL.GL_LIGHTING);  // shouldn't need this
                if (surfDisplayListNum > 0)
                {
                    gl.glDeleteLists(surfDisplayListNum, 1);
                    surfDisplayListNum = 0;
                }
                drawSurface(glPanel3d, displayLayers, displayProperties, propertyVolume, propertyChanged);
            }
            propertyChanged = false;
            colorChanged = false;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "displaySurface", e);
        }
        finally
        {
            // glPanel3d.resetViewShift(gl);
            gl.glFrontFace(GL.GL_CCW);  // restore default
        }
    }

    private void drawSurface(StsGLPanel3d glPanel3d, boolean displayLayers, boolean displayProperties, StsPropertyVolume propertyVolume, boolean propertyChanged)
    {
        GL gl = glPanel3d.getGL();
        if (gl == null) return;
        gl.glEnable(GL.GL_LIGHTING); // Shouldn't need this
        GLU glu = glPanel3d.getGLU();

        if (grid instanceof StsXYSurfaceGridable)
        {
            StsXYSurfaceGridable xyGrid = (StsXYSurfaceGridable) grid;
            float[][] pointsZ = xyGrid.getPointsZ();
            drawXYQuadStrips(gl, glu, xyGrid, pointsZ);
        }
        else if(!displayLayers)
            drawQuadStrips(gl, glu);
        else
            drawQuadStripPolygons(gl, glu, displayLayers, displayProperties, propertyVolume);

//		gl.glFrontFace(GL.GL_CW);

        drawPolygons(gl, glu, displayLayers, displayProperties, propertyVolume, propertyChanged);
    }

    private void displayContours(StsModel model, StsGLPanel3d glPanel3d, boolean displayGaps, boolean displaySubPolygons, boolean displayFill, StsColor color, StsPropertyVolume propertyVolume)
    {
        GL gl = glPanel3d.getGL();
        if (gl == null) return;
        try
        {
            StsColor gridColor = displayFill ? StsColor.BLACK : color;
            gridColor.setGLColor(gl);

            boolean useDisplayLists = model.useDisplayLists;

            if (contourDisplayListNum == 0 && useDisplayLists)  // build display list
            {
                contourDisplayListNum = gl.glGenLists(1);
                if (contourDisplayListNum == 0)
                {
                    StsMessageFiles.logMessage("System Error in StsGrid.displaySurface: " +
                        "Failed to allocate a display list");
                    return;
                }

                gl.glNewList(contourDisplayListNum, GL.GL_COMPILE_AND_EXECUTE);
                drawContours(glPanel3d, displayGaps, displaySubPolygons, propertyVolume);
                //            normals = null; // delete normals after building for display list
                gl.glEndList();

                //timer.stop("display list surface setup: ");

            }
            else if (useDisplayLists) // use existing display list
            {
                // gl.glEnable(GL.GL_LIGHTING);  // shouldn't need this
                gl.glCallList(contourDisplayListNum);
                //timer.stop("display list surface draw: ");
            }
            else // immediate mode draw
            {
                gl.glEnable(GL.GL_LIGHTING);  // shouldn't need this
                if (contourDisplayListNum > 0)
                {
                    gl.glDeleteLists(contourDisplayListNum, 1);
                    contourDisplayListNum = 0;
                }
                drawContours(glPanel3d, displayGaps, displaySubPolygons, propertyVolume);
            }

        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "displayContours", e);
        }
        finally
        {
            // glPanel3d.resetViewShift(gl);
            gl.glFrontFace(GL.GL_CCW);  // restore default
        }
    }

    private void drawContours(StsGLPanel3d glPanel3d, boolean displaySubPolygons, boolean displayProperties, StsPropertyVolume propertyVolume)
    {
        GL gl = glPanel3d.getGL();
        if (gl == null) return;
        gl.glEnable(GL.GL_LIGHTING); // Shouldn't need this
        GLU glu = glPanel3d.getGLU();

        if (displaySubPolygons)
        {
            drawQuadStripContours(gl, glu, displayProperties, propertyVolume);
        }
        else if (grid instanceof StsXYSurfaceGridable)
        {
            StsXYSurfaceGridable xyGrid = (StsXYSurfaceGridable) grid;
            float[][] pointsZ = xyGrid.getPointsZ();
            drawXYQuadStripContours(gl, glu, xyGrid, pointsZ);
        }
        else
            drawQuadStripContours(gl, glu);

//		gl.glFrontFace(GL.GL_CW);

        drawPolygonContours(gl, glu, displaySubPolygons);
    }

    private void setViewShift(StsGLPanel3d glPanel3d, GL gl)
    {
        if (grid instanceof StsSection || grid instanceof StsCursorSection)
            glPanel3d.setViewShift(gl, StsGraphicParameters.sectionShift);
    }

    private void drawQuadStripPolygons(GL gl, GLU glu, boolean displaySubPolygons, boolean displayProperties, StsPropertyVolume propertyVolume)
    {
        QuadStrip q = null;

        try
        {
            if (quadStrips == null) return;
            String name = toString();
            int nStrips = quadStrips.getSize();
            for (int n = 0; n < nStrips; n++)
            {
                q = (QuadStrip) quadStrips.getElement(n);
                if (q.polygons != null)
                    drawQuadStripPolygons(q, gl, glu, displaySubPolygons, displayProperties, propertyVolume, propertyChanged);
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.drawQuadStrips() failed." +
                q.rowNumber + q.firstCol, e, StsException.WARNING);
        }
    }

    private void drawQuadStripContours(GL gl, GLU glu, boolean displayProperties, StsPropertyVolume propertyVolume)
    {
        QuadStrip q = null;
        float[] point, normal;
        int i = -1, j = -1;

        try
        {
            if (quadStrips == null) return;
            String name = toString();
            int nStrips = quadStrips.getSize();
            for (int n = 0; n < nStrips; n++)
            {
                q = (QuadStrip) quadStrips.getElement(n);
                if (q.polygons != null)
                    drawQuadStripContours(q, gl, glu, displayProperties, propertyVolume, propertyChanged);
            }
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "drawQuadStripContours" + q.rowNumber + q.firstCol, e);
        }
    }

    private void drawXYQuadStrips(GL gl, GLU glu, StsXYSurfaceGridable xyGrid, float[][] pointsZ)
    {
        QuadStrip q = null;
        float[] point, normal;
        int i = -1, j = -1;

        try
        {
            if (quadStrips == null) return;
            float xMin = xyGrid.getXMin();
            float yMin = xyGrid.getYMin();
            float xInc = xyGrid.getXInc();
            float yInc = xyGrid.getYInc();
            int rowMin = xyGrid.getRowMin();
            int colMin = xyGrid.getColMin();
            int nStrips = quadStrips.getSize();
            for (int n = 0; n < nStrips; n++)
            {
                q = (QuadStrip) quadStrips.getElement(n);
                i = q.rowNumber;
                j = q.firstCol;
                float x = xMin + j * xInc;
                float y = yMin + i * yInc;
                float y1 = y + yInc;

                gl.glBegin(GL.GL_QUAD_STRIP);
                for (j = q.firstCol; j <= q.lastCol; j++)
                {
                    float[] gridNormal = grid.getNormal(i, j);
                    /*
                         if(gridNormal == null)
                         {
                             StsException.systemError("Bad normal at row " + i + " col " + j);
                             continue;

                         }
                     */
                    gl.glNormal3fv(grid.getNormal(i, j), 0);
                    gl.glVertex3f(x, y, pointsZ[i - rowMin][j - colMin]);
                    gl.glNormal3fv(grid.getNormal(i + 1, j), 0);
                    gl.glVertex3f(x, y1, pointsZ[i + 1 - rowMin][j - colMin]);
                    x += xInc;
                }
                gl.glEnd();
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.drawQuadStrips() failed." +
                " grid " + xyGrid.toString() + " row: " + i + " col: " + j, e, StsException.WARNING);
        }
    }

    private void drawQuadStrips(GL gl, GLU glu)
    {
        QuadStrip q = null;
        float[] point, normal;
        int i = -1, j = -1;

        try
        {
            if (quadStrips == null) return;

            int nStrips = quadStrips.getSize();
            for (int n = 0; n < nStrips; n++)
            {
                q = (QuadStrip) quadStrips.getElement(n);
                i = q.rowNumber;
                gl.glBegin(GL.GL_QUAD_STRIP);
                for (j = q.firstCol; j <= q.lastCol; j++)
                {
                    gl.glNormal3fv(grid.getNormal(i, j), 0);
                    gl.glVertex3fv(grid.getXYZorT(i, j), 0);
                    gl.glNormal3fv(grid.getNormal(i + 1, j), 0);
                    gl.glVertex3fv(grid.getXYZorT(i + 1, j), 0);
                }
                gl.glEnd();
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.drawQuadStrips() failed." +
                "row: " + i + " col: " + j, e, StsException.WARNING);
        }
    }

    private void drawQuadStripPolygons(QuadStrip q, GL gl, GLU glu, boolean displaySubPolygons, boolean displayProperties, StsPropertyVolume propertyVolume, boolean propertyChanged)
    {
        StsException stsException = null;
        StsPolygon polygon = null;
        int n = -1;
        try
        {
            if (q.polygons == null) return;

            int firstCol = q.firstCol;
            int lastCol = q.lastCol;
            int nCols = lastCol - firstCol;
            for (int j = 0; j < nCols; j++)
            {
                Object[] cellPolygons = q.polygons[j];
                if (cellPolygons == null) continue;
                int nPolygons = cellPolygons.length;
                for (n = 0; n < nPolygons; n++)
                {
                    try
                    {
                        polygon = (StsPolygon) cellPolygons[n];
                        if (polygon != null) polygon.draw(gl, glu, isPlanar, normal, displaySubPolygons, displayProperties, propertyVolume, this.propertyChanged);
                    }
                    catch (Exception e)
                    {
                        if (stsException == null) stsException = new StsException(StsException.WARNING,
                            "StsEdgeLoop.drawQuadStripPolygons() failed.",
                            q.toString() + " polygon " + n + ": " + polygon.getLabel());
                    }
                }
            }
        }
        catch (Exception e)
        {
            if (stsException == null)
            {
                String polygonLabel = "NULL";
                if (polygon != null) polygonLabel = polygon.getLabel();
                stsException = new StsException(StsException.WARNING,
                    "StsEdgeLoop.drawQuadStripPolygons() failed.",
                    q.toString() + " polygon " + n + ": " + polygonLabel);
            }
        }
        finally
        {
            if (stsException != null) StsException.outputException(stsException);
        }
    }
   private void drawQuadStripContours(QuadStrip q, GL gl, GLU glu, boolean displayProperties, StsPropertyVolume propertyVolume, boolean propertyChanged)
    {
        StsException stsException = null;
        StsPolygon polygon = null;
        int n = -1;
        try
        {
            if (q.polygons == null) return;

            int firstCol = q.firstCol;
            int lastCol = q.lastCol;
            int nCols = lastCol - firstCol;
            for (int j = 0; j < nCols; j++)
            {
                Object[] cellPolygons = q.polygons[j];
                if (cellPolygons == null) continue;
                int nPolygons = cellPolygons.length;
                for (n = 0; n < nPolygons; n++)
                {
                    try
                    {
                        polygon = (StsPolygon) cellPolygons[n];
                        if (polygon != null) polygon.draw(gl, glu, isPlanar, normal, true, displayProperties, propertyVolume, propertyChanged);
                    }
                    catch (Exception e)
                    {
                        if (stsException == null) stsException = new StsException(StsException.WARNING,
                            "StsEdgeLoop.drawQuadStripPolygons() failed.",
                            q.toString() + " polygon " + n + ": " + polygon.getLabel());
                    }
                }
            }
        }
        catch (Exception e)
        {
            if (stsException == null)
            {
                String polygonLabel = "NULL";
                if (polygon != null) polygonLabel = polygon.getLabel();
                stsException = new StsException(StsException.WARNING,
                    "StsEdgeLoop.drawQuadStripPolygons() failed.",
                    q.toString() + " polygon " + n + ": " + polygonLabel);
            }
        }
        finally
        {
            if (stsException != null) StsException.outputException(stsException);
        }
    }

    private void drawXYQuadStripContours(GL gl, GLU glu, StsXYSurfaceGridable xyGrid, float[][] pointsZ)
    {
        QuadStrip q = null;
        float[] point, normal;
        int i = -1, j = -1;

        try
        {
            if (quadStrips == null) return;
            float xMin = xyGrid.getXMin();
            float yMin = xyGrid.getYMin();
            float xInc = xyGrid.getXInc();
            float yInc = xyGrid.getYInc();
            int rowMin = xyGrid.getRowMin();
            int colMin = xyGrid.getColMin();
            int nStrips = quadStrips.getSize();
            for (int n = 0; n < nStrips; n++)
            {
                q = (QuadStrip) quadStrips.getElement(n);
                i = q.rowNumber;
                j = q.firstCol;
                float x = xMin + j * xInc;
                float y = yMin + i * yInc;
                float y1 = y + yInc;

                gl.glBegin(GL.GL_QUAD_STRIP);
                for (j = q.firstCol; j <= q.lastCol; j++)
                {
                    float[] gridNormal = grid.getNormal(i, j);
                    /*
                         if(gridNormal == null)
                         {
                             StsException.systemError("Bad normal at row " + i + " col " + j);
                             continue;

                         }
                     */
                    gl.glNormal3fv(grid.getNormal(i, j), 0);
                    gl.glVertex3f(x, y, pointsZ[i - rowMin][j - colMin]);
                    gl.glNormal3fv(grid.getNormal(i + 1, j), 0);
                    gl.glVertex3f(x, y1, pointsZ[i + 1 - rowMin][j - colMin]);
                    x += xInc;
                }
                gl.glEnd();
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.drawQuadStrips() failed." +
                " grid " + xyGrid.toString() + " row: " + i + " col: " + j, e, StsException.WARNING);
        }
    }

    private void drawQuadStripContours(GL gl, GLU glu)
    {
        QuadStrip q = null;
        StsPoint[] points = new StsPoint[4];
        int i = -1, j = -1;

        try
        {
            if (quadStrips == null) return;

            int nStrips = quadStrips.getSize();
            for (int n = 0; n < nStrips; n++)
            {
                q = (QuadStrip) quadStrips.getElement(n);
                i = q.rowNumber;
                for (j = q.firstCol; j <= q.lastCol; j++)
                {
                    points[0] = grid.getPoint(i, j);
                    points[1] = grid.getPoint(i+1, j);
                    points[2] = grid.getPoint(i+1, j+1);
                    points[3] = grid.getPoint(i, j+1);
                }
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsEdgeLoop.drawQuadStrips() failed." +
                "row: " + i + " col: " + j, e, StsException.WARNING);
        }
    }

    private void drawPolygonContours(GL gl, GLU glu, boolean drawSubPolygons)
    {
        if (polygons == null) return;

//        StsPolygon.setLogFile(model.getLogFile());

        int nPolygons = polygons.getSize();
        String name = grid.toString();
        for (int n = 0; n < nPolygons; n++)
        {
            StsPolygon polygon = (StsPolygon) polygons.getElement(n);
            if (polygon != null)
            {
                if (!polygon.drawContours(gl, glu, isPlanar, normal, drawSubPolygons)) polygon = null;
            }
        }
    }

    private void drawPolygons(GL gl, GLU glu, boolean displaySubPolygons, boolean displayProperties, StsPropertyVolume propertyVolume, boolean propertyChanged)
    {
        if (polygons == null) return;

//        StsPolygon.setLogFile(model.getLogFile());

        int nPolygons = polygons.getSize();
        String name = grid.toString();
        for (int n = 0; n < nPolygons; n++)
        {
            StsPolygon polygon = (StsPolygon) polygons.getElement(n);
            if (polygon != null)
            {
                if (!polygon.draw(gl, glu, isPlanar, normal, displaySubPolygons, displayProperties, propertyVolume, propertyChanged)) polygon = null;
            }
        }
    }

    private void displayGrid(StsModel model, StsGLPanel3d glPanel3d, boolean displayFill, StsColor color)
    {
        GL gl = glPanel3d.getGL();
        if (gl == null) return;
        gl.glDisable(GL.GL_LIGHTING);

        StsColor gridColor = displayFill ? StsColor.BLACK : color;
        gridColor.setGLColor(gl);

        boolean useDisplayLists = model.getBooleanProperty("Use Display Lists");
        if (StsGLPanel.debugProjectionMatrix)
            StsException.systemDebug("viewShift called for displayGrid of " + toString());
        glPanel3d.setViewShift(gl, StsGraphicParameters.gridShift);

        if (gridDisplayListNum == 0 && useDisplayLists)  // build display list
        {
            gridDisplayListNum = gl.glGenLists(1);
            if (gridDisplayListNum == 0)
            {
                StsMessageFiles.logMessage("System Error in StsEdgeLoop.displayGrid: " +
                    "Failed to allocate a display list");
                return;
            }

            gl.glNewList(gridDisplayListNum, GL.GL_COMPILE_AND_EXECUTE);
            drawGridLines(gl);
            gl.glEndList();

            //timer.stop("display list surface setup: ");

        }
        else if (useDisplayLists) // use existing display list
        {
            gl.glCallList(gridDisplayListNum);
            //timer.stop("display list surface draw: ");
        }
        else
        {
            if (gridDisplayListNum > 0)
            {
                gl.glDeleteLists(gridDisplayListNum, 1);
                gridDisplayListNum = 0;
            }
            drawGridLines(gl);
        }
        gl.glEnable(GL.GL_LIGHTING);
        glPanel3d.resetViewShift(gl);
    }

    private void drawGridLines(GL gl)
    {
        if (rowGridLines != null) drawGridLines(rowGridLines, gl);
        if (colGridLines != null) drawGridLines(colGridLines, gl);
    }

    // draw all gridLines for this array of gridLines (rowGridLines or colGridLines)
    private void drawGridLines(GridArray2D gridLines, GL gl)
    {
        StsEdgeLinkable gridLine;

        try
        {
            if (gridLines == null) return;

            gl.glLineWidth(StsGraphicParameters.gridLineWidth);
            Iterator iter = gridLines.iterator();
            while (iter.hasNext())
            {
                gridLine = (StsEdgeLinkable) iter.next();
                drawGridLine(gridLine, gl);
            }
        }
        catch (Exception e){} // mainDebug, so don't process exception
    }

    // get gridLines for this specific row or col
    public Object[] getGridLines(int rowOrCol, int rowCol)
    {
        if (rowOrCol == ROW && rowGridLines != null)
            return rowGridLines.getList(rowCol);
        else if (rowOrCol == COL && colGridLines != null)
            return colGridLines.getList(rowCol);
        else
            return null;
    }

    // draw gridLines for this specific row or col
    public void drawGridLines(int rowOrCol, int rowCol, GL gl)
    {
        Object[] gridLines = getGridLines(rowOrCol, rowCol);
        if (gridLines == null) return;
        for (int n = 0; n < gridLines.length; n++)
        {
            StsEdgeLinkable gridLine = (StsEdgeLinkable) gridLines[n];
            drawGridLine(gridLine, gl);
        }
    }

    private void drawGridLine(StsEdgeLinkable gridLine, GL gl)
    {
        float[][] xyzPoints = gridLine.getXYZPoints();
        StsGLDraw.drawSegmentedLine(gl, xyzPoints);
    }

    /*
        private void drawGridLines(GL gl)
        {
            StsGridLine gridLine;

            try
            {
                if(gridLines == null) return;

                gl.glLineWidth(StsGraphicParameters.gridLineWidth);
                int nGridLines = gridLines.getSize();
                for(int n = 0; n < nGridLines; n++)
                {
                    gridLine = (StsGridLine)gridLines.getElement(n);
                    drawGridLine(gridLine, gl);
                }
            }
            catch(Exception e) {} // mainDebug, so don't process exception
        }

        private void drawGridLine(StsGridLine gridLine, GL gl)
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
                            xyz = grid.getXYZ(row, col);
                            if(xyz != null)  gl.glVertex3fv(xyz);
                        }
                    }
                    else
                    {
                        col = gridLine.getRowCol();
                        for(row = first; row <= last; row++)
                        {
                            xyz = grid.getXYZ(row, col);
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

                StsException.outputException("StsEdgeLoop.debugDrawGridLine() failed. " +
                    StsParameters.rowCol(rowOrCol) + rowCol, e, StsException.WARNING);
            }
        }
    */
    public void deleteDisplayLists(GL gl)
    {
        deleteSurfDisplayList(gl);
        deleteGridDisplayList(gl);
    }

    public void deleteSurfDisplayList(GL gl)
    {
        if (surfDisplayListNum <= 0) return;
        gl.glDeleteLists(surfDisplayListNum, 1);
        surfDisplayListNum = 0;
    }

    public void deleteGridDisplayList(GL gl)
    {
        if (gridDisplayListNum <= 0) return;
        gl.glDeleteLists(gridDisplayListNum, 1);
        gridDisplayListNum = 0;
    }

    public void deleteTransientArrays()
    {
        loopLinks = null;
        linkedGrid = null;
        cellTypeGrid = null;
    }

    public void checkSetPropertyChanged(String newDisplayMode)
    {
        if(newDisplayMode == displayMode) return;
        colorChanged = true;
        if(newDisplayMode == StsBuiltModelClass.displayLayersString || newDisplayMode == StsBuiltModelClass.staticGetDisplayModeString())
            propertyChanged = true;
        displayMode = newDisplayMode;
    }

    public void zoneSidePropertyChanged(StsPropertyVolume propertyVolume)
    {
        propertyChanged = true;
        Iterator gridPolygonIterator = getZoneSideLayerPolygonIterator();
        if (gridPolygonIterator == null) return;
        while (gridPolygonIterator.hasNext())
        {
            StsPolygon polygon = (StsPolygon) gridPolygonIterator.next();
            polygon.setPropertyColor(propertyVolume);
        }
    }

    public void propertyChanged(StsPropertyVolume propertyVolume)
    {
        propertyChanged = true;
        Iterator gridPolygonIterator = new GridPolygonIterator();
        if (gridPolygonIterator == null) return;
        while (gridPolygonIterator.hasNext())
        {
            StsPolygon polygon = (StsPolygon) gridPolygonIterator.next();
            polygon.setPropertyColor(propertyVolume);
        }
    }

    public Iterator getZoneSideGridPolygonIterator()
    {
        if (polygonType == StsPolygon.GRID)
            return polygons.iterator();
        else if (polygonType == StsPolygon.SECTION)
            return new GridPolygonIterator();
        else
            return StsList.emptyIterator();
    }

    public Iterator getZoneSideLayerPolygonIterator()
    {
        if (polygonType == StsPolygon.GRID)
            return new LayerPolygonIterator();
        else if (polygonType == StsPolygon.SECTION)
            return new GridPolygonIterator();
        else
            return StsList.emptyIterator();
    }

    class LayerPolygonIterator implements Iterator
    {
        Iterator<StsPolygon> gridIterator = null;
        Iterator layerIterator = StsPolygon.emptyIterator();

        LayerPolygonIterator()
        {
            gridIterator = getPolygonIterator();
        }

        public boolean hasNext()
        {
            if(layerIterator.hasNext()) return true;
            if (gridIterator.hasNext())
            {
                StsPolygon polygon = gridIterator.next();
                layerIterator = polygon.getSubPolygonIterator();
                return layerIterator.hasNext();
            }
            return false;
        }

        public Object next()
        {
            return layerIterator.next();
        }

        public void remove()
        {
            StsException.systemError(this, "remove", "Not allowed to remove a quad/polygon element from iterator.");
        }
    }
    /**
     * This iterator is for a non-vertical section.  The edgeLoop is subdivided into quad cells with trim polygons on the edges.
     * We iterate over each quad and each trim polygons.  Each of these in turn will contain GRID subPolygons.
     */
    class GridPolygonIterator implements Iterator<StsPolygon>
    {
        Iterator<StsPolygon> iterator = null;
        boolean polygonsStarted = false;

        GridPolygonIterator()
        {
            if (quadStrips != null)
            {
                iterator = new QuadCellGridPolygonIterator();
            }
            else if (polygons != null)
            {
                iterator = new EdgePolygonGridPolygonIterator();
                polygonsStarted = true;
            }
            else
                iterator = StsPolygon.emptyIterator();
        }

        public boolean hasNext()
        {
            // return current iterator status
            if (iterator.hasNext())
                return true;
            // if this iterator is finished and polygons aren't started, make it the iterator
            if (!polygonsStarted && polygons != null)
            {
                iterator = new EdgePolygonGridPolygonIterator();
                polygonsStarted = true;
                return iterator.hasNext();
            }
            // iteration over quads and polygons is completed
            return false;
        }

        public StsPolygon next()
        {
            return iterator.next();
        }

        public void remove()
        {
            StsException.systemError(this, "remove", "Not allowed to remove a quad/polygon element from iterator.");
        }
    }

    class QuadCellLayerPolygonIterator implements Iterator<StsPolygon>
    {
        Iterator<QuadStrip> quadStripIterator = quadStrips.iterator();
        Iterator<StsPolygon> polygonIterator = StsPolygon.emptyIterator();
        Iterator<StsPolygon> layerIterator = StsPolygon.emptyIterator();

        QuadCellLayerPolygonIterator()
        {
            if(quadStrips != null)
                quadStripIterator = quadStrips.iterator();
            else
                quadStripIterator = QuadStrip.emptyIterator();
        }

        public boolean hasNext()
        {
            if(layerIterator.hasNext())
                return true;
            if (polygonIterator.hasNext())
            {
                StsPolygon polygon = polygonIterator.next();
                layerIterator = polygon.getSubPolygonIterator();
                return layerIterator.hasNext();
            }
            if (quadStripIterator.hasNext())
            {
                QuadStrip quadStrip = quadStripIterator.next();
                polygonIterator = quadStrip.getPolygonIterator();
                if (polygonIterator.hasNext())
                {
                    StsPolygon polygon = polygonIterator.next();
                    layerIterator = polygon.getSubPolygonIterator();
                    return layerIterator.hasNext();
                }
                else
                    return false;
            }
            else
                return false;
        }

        public StsPolygon next()
        {
            return polygonIterator.next();
        }

        public void remove(){ StsException.systemError(this, "remove", "Cannot remove element from iterator."); }
    }

    /** This iterator iterates over each quadCell and over the polygons within these cells which are GRID polygons. */
    class QuadCellGridPolygonIterator implements Iterator<StsPolygon>
    {
        Iterator<QuadStrip> quadStripIterator;
        Iterator<StsPolygon> polygonIterator = StsPolygon.emptyIterator();

        QuadCellGridPolygonIterator()
        {
           if(quadStrips != null)
                quadStripIterator = quadStrips.iterator();
            else
                quadStripIterator = QuadStrip.emptyIterator();
        }

        public boolean hasNext()
        {
            if (polygonIterator.hasNext())
                return true;
            if (quadStripIterator.hasNext())
            {
                QuadStrip quadStrip = quadStripIterator.next();
                polygonIterator = quadStrip.getPolygonIterator();
                return polygonIterator.hasNext();
            }
            else
                return false;
        }

        public StsPolygon next()
        {
            return polygonIterator.next();
        }

        public void remove(){ StsException.systemError(this, "remove", "Cannot remove element from iterator."); }
    }

    class EdgePolygonGridPolygonIterator implements Iterator<StsPolygon>
    {
        Iterator<StsPolygon> edgePolygonIterator = StsPolygon.emptyIterator();
        Iterator<StsPolygon> gridPolygonIterator = StsPolygon.emptyIterator();
        boolean polygonsStarted = false;

        EdgePolygonGridPolygonIterator()
        {
            if (polygons != null)
                edgePolygonIterator = polygons.iterator();
            else
                edgePolygonIterator = StsPolygon.emptyIterator();
        }

        public boolean hasNext()
        {
            if (gridPolygonIterator.hasNext())
                return true;
            if (edgePolygonIterator.hasNext())
            {
                StsPolygon edgePolygon = edgePolygonIterator.next();
                gridPolygonIterator = edgePolygon.getSubPolygonIterator();
                return hasNext();
            }
            else
                return false;
        }

        public StsPolygon next()
        {
            return gridPolygonIterator.next();
        }

        public void remove(){ StsException.systemError(this, "remove", "Cannot remove element from iterator."); }
    }

    class EdgePolygonLayerPolygonIterator implements Iterator<StsPolygon>
    {
        Iterator<StsPolygon> edgePolygonIterator = StsPolygon.emptyIterator();
        Iterator<StsPolygon> polygonIterator = StsPolygon.emptyIterator();
        Iterator<StsPolygon> layerIterator = StsPolygon.emptyIterator();

        boolean polygonsStarted = false;

        EdgePolygonLayerPolygonIterator()
        {
            if (polygons != null)
                edgePolygonIterator = polygons.iterator();
        }

        public boolean hasNext()
        {
            if(layerIterator.hasNext())
                return true;
            if (polygonIterator.hasNext())
            {
                StsPolygon polygon = polygonIterator.next();
                layerIterator = polygon.getSubPolygonIterator();
                return layerIterator.hasNext();
            }
            if (edgePolygonIterator.hasNext())
            {
                StsPolygon edgePolygon = edgePolygonIterator.next();
                polygonIterator = edgePolygon.getSubPolygonIterator();
                if (polygonIterator.hasNext())
                {
                    StsPolygon polygon = polygonIterator.next();
                    layerIterator = polygon.getSubPolygonIterator();
                    return layerIterator.hasNext();
                }
                else
                    return false;
            }
            else
                return false;
        }

        public StsPolygon next()
        {
            return polygonIterator.next();
        }

        public void remove(){ StsException.systemError(this, "remove", "Cannot remove element from iterator."); }
    }

    public Iterator getPolygonIterator()
    {
        ArrayList<StsPolygon> polygonList = new ArrayList<StsPolygon>();
        if (polygons != null)
        {
            int nPolygons = polygons.getSize();
            for(int n = 0; n < nPolygons; n++)
                polygonList.add((StsPolygon)polygons.getElement(n));

        }
        return polygonList.iterator();
    }

    public StsGridBoundingBox getGridBoundingBox()
    {
        if(linkedGrid == null) return null;
        return linkedGrid.getLoopBoundingBox();
    }

    /*
        private void drawQuadStripProps(GL gl, StsPropertyVolumeOld pv, boolean displayGaps)
        {
            QuadStrip q = null;
            float[] p00, p01, p10, p11, n00, n01, n10, n11;
            int i = -1, j = -1, ii, jj;

            try
            {
                gridRect = new StsGridRectangle(this);
                int propLayerNumber = surface.getPropertyLayerNumber();
                StsColor colors[][] = pv.getSliceColors(gridRect, propLayerNumber);

                int nStrips = quadStrips.getSize();
                for (int n = 0; n < nStrips; n++)
                {
                    q = (QuadStrip)quadStrips.getElement(n);
                    if(!displayGaps && q.cellType == CELL_GAP) continue;

                    i = q.rowOrColNumber;
                    ii = i - pointRowMin;
                    j = q.firstCol;
                    jj = j - pointColMin;

                    p01 = getXYZ(i, j);
                    p11 = getXYZ(i+1, j);
                    n01 = getNormal(i, j);
                    n11 = getNormal(i+1, j);

                    gl.glBegin(GL.GL_QUADS);
                    for(; j < q.lastCol; j++, jj++)
                    {
                        n00 = n01;
                        n10 = n11;
                        p00 = p01;
                        p10 = p11;

                        StsColor.setGLColor(gl, colors[ii][jj]);

                        if(n00 != null) gl.glNormal3fv(n00);
                        if(p00 != null)  gl.glVertex3fv(p00);
                        if(n10 != null)  gl.glNormal3fv(n10);
                        if(p10 != null)  gl.glVertex3fv(p10);

                        n01 = getNormal(i, j+1);
                        p01 = getXYZ(i, j+1);
                        n11 = getNormal(i+1, j+1);
                        p11 = getXYZ(i+1, j+1);

                        if(n11 != null)  gl.glNormal3fv(n11);
                        if(p11 != null)  gl.glVertex3fv(p11);
                        if(n01 != null) gl.glNormal3fv(n01);
                        if(p01 != null)  gl.glVertex3fv(p01);
                    }
                    gl.glEnd();
                }
            }
            catch(Exception e)
            {
                StsException.outputException("StsEdgeLoop.drawSurface1() failed." +
                    "row: " + i + " col: " + j, e, StsException.WARNING);
            }
        }

        private void drawPolygonProps(GL gl, boolean displayGaps)
        {
            GLU glu = currentModel.win3d.glPanel3d.getGLU();

            int nPolygons = polygons.getSize();
            for(int n = 0; n < nPolygons; n++)
            {
                StsPolygon polygon = (StsPolygon)polygons.getElement(n);
                StsColor color = surface.getPropColor(polygon.getRow(), polygon.getCol());
                if(displayGaps || polygon.getType() != FAULT) polygon.draw(gl, glu, this, color);
            }
        }

        private void checkQuadStripNormals()
        {
            QuadStrip q = null;
            int i = -1, j = -1;

            try
            {
                if(quadStrips == null) return;

                int nStrips = quadStrips.getSize();
                for (int n = 0; n < nStrips; n++)
                {
                    q = (QuadStrip)quadStrips.getElement(n);
                    i = q.rowOrColNumber;
                    for(j = q.firstCol; j <= q.lastCol; j++)
                    {
                        checkNormal(i, j);
                        checkNormal(i+1, j);
                    }
                }
            }
            catch(Exception e)
            {
                StsException.outputException("StsEdgeLoop.checkQuadStripNormals() failed." +
                    "row: " + i + " col: " + j, e, StsException.WARNING);
            }
        }

        public final void checkPointZ(int row, int col)
        {
            if(pointsZ[row-pointRowMin][col-pointColMin] != nullValue) return;

            StsException.systemError(getLabel() +
                " filling in a pointZ at row: " + row + " col: " + col);
            pointsZ[row-pointRowMin][col-pointColMin] = surfaceGrid.getPointZ(row, col);
        }

        public final void checkNormal(int row, int col)
        {
            if(normals[row-pointRowMin][col-pointColMin] != null) return;

            StsException.systemError(getLabel() + " at row: " + row +
                " col: " + col + " needs a normal.");
            normals[row-pointRowMin][col-pointColMin] = vertNormal;
        }

    */
    public void debugDisplay(StsGLPanel3d glPanel3d)
    {
        if (loopLinks == null) return;

        GL gl = glPanel3d.getGL();
        if (gl == null) return;

        gl.glDisable(GL.GL_LIGHTING);
        glPanel3d.setViewShift(gl, StsGraphicParameters.vertexShift);

        int nLinks = loopLinks.getSize();
        for (int n = 0; n < nLinks; n++)
        {
            StsEdgeLoopRadialGridLink link = (StsEdgeLoopRadialGridLink) loopLinks.getElement(n);
            link.debugDisplay(glPanel3d, 0.0f);
        }
        gl.glEnable(GL.GL_LIGHTING);
        glPanel3d.resetViewShift(gl);
    }


    public void propertyTypeChanged()
    {
        propertyChanged = true;
    }

    public String getLabel()
    {
        if (grid == null) return new String("null");
        else return grid.getLabel();
    }

    public String toString()
    {
        if (grid == null) return super.toString();
        else return grid.toString();
    }
}

