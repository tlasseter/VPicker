
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Build;

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

import java.util.*;

public class StsIntersectSectionEdges extends StsAction implements Runnable
{
    private boolean actionOK = false;

    StsRotatedGridBoundingBox grid;
    int nRows, nCols;
    private EdgeCell[][] edgeCells;
    ArrayList intersections;
    StsSectionEdge[] differentEdges = null; // scratch: used to check every cell, so allocate once and resize as needed

    HashMap intersectedEdges; // a set of EdgeIntersections with the associated StsSectionEdge as the hashCode

    static public final float maxGridGap = 2.0f;

    // type of intersection
    public static final int MINUS = StsParameters.MINUS;
    public static final int PLUS = StsParameters.PLUS;
    public static final int NONE = StsParameters.NONE;
    public static final int CROSS = StsParameters.PLUS_AND_MINUS;

   	/** Makes intersections between sectionEdges: join if they are in the same cell,
        splitting/joining if they are in same cell and all ends are longer than two cells:
        shorter end if two cross */
    public StsIntersectSectionEdges(StsActionManager actionManager)
    {
        super(actionManager, true);
    }

	public void run()
    {

		if(!actionOK()) return;

        try
        {
            model.disableDisplay();
            statusArea.setTitle("Intersecting section edges: ");
            intersectEdges();
            actionManager.endCurrentAction();

            model.enableDisplay();
            model.win3d.win3dDisplay();
//            model.glPanel.getActionManager().fireChangeEvent();
        }
        catch(Exception e)
        {
            actionOK = false;
            return;
        }
    }

	private boolean actionOK()
    {
    	try
        {
            actionOK = true;

            if (!model.hasObjectsOfType(StsModelSurface.class, StsModelSurface.MODEL))
            {
                logMessage("No horizons have been built." +  "  Terminating action.");
                actionOK = false;
            }
            else
            {
                StsClass sectionEdges = model.getCreateStsClass(StsSectionEdge.class);
                if(sectionEdges.getSize() <= 0)
                {
                    logMessage("There are no section edges currently."
                        +  "  Terminating action.");
                    actionOK = false;
                }
            }
	        return actionOK;
        }
        catch(Exception e)
        {
            actionOK = false;
            return actionOK;
        }
    }

    public boolean end()
    {
        if (!actionOK) logMessage("Fault cuts from imported cuts failed.");
        return actionOK;
    }

    boolean insideGrid(int row, int col)
    {
        return row >= 0 && row < nRows && col >= 0 && col < nCols;
    }

    private void intersectEdges()
    {
        try
        {
            constructGrid();
            edgeCells = new EdgeCell[nRows][nCols];

            StsObject[] sections = model.getObjectList(StsSection.class);
            int nSectionEdges = sections.length;
            intersectedEdges = new HashMap(nSectionEdges);

            for(int n = 0; n < nSectionEdges; n++)
            {
                StsSection section = (StsSection)sections[n];
                StsSectionEdge edge = (StsSectionEdge)section.getSectionEdges().getFirst();
                addEdge(edge);
            }

            createIntersectionList();
            makeIntersections();
        }
        catch(Exception e)
        {
            StsException.outputException("StsIntersectSectionEdges.intersectEdges() failed.",
                e, StsException.WARNING);
        }
    }

    static int maxNCells = 20;
    private void constructGrid()
    {
        StsRotatedGridBoundingBox projectGrid = model.getProject().getRotatedBoundingBox();
        grid = new StsRotatedGridBoundingBox(projectGrid, false);
        nRows = grid.getNRows();
        int maxNRows = maxNCells + 1;
        if(nRows > maxNRows)
        {
            grid.setNRows(maxNRows);
            nRows = maxNRows;
            float yInc = grid.getYSize()/(maxNCells);
            grid.setYInc(yInc);
        }
        nCols = grid.getNCols();
        int maxNCols = maxNCells + 1;
        if(nCols > maxNCols)
        {
            grid.setNCols(maxNCols);
            nCols = maxNCols;
            float xInc = grid.getXSize()/maxNCells;
            grid.setXInc(xInc);
        }
    }

    private void addEdge(StsSectionEdge edge)
    {
        StsGridCrossingPoint point0, point1, gridPoint0, gridPoint1;
        float rowF0, colF0, rowF1, colF1;
        int row, col, p, g;

        try
        {
            edge.setPointIndexFs();
            IntersectedEdge intersectedEdge = new IntersectedEdge(edge);
            Integer key = new Integer(edge.getIndex());
            intersectedEdges.put(key, intersectedEdge);

            StsPoint[] points = edge.getPoints();
            int edgeIndex = edge.getIndex();
            int nPoints = points.length;

            ArrayList cellSegments = new ArrayList(100);
            point1 = new StsGridCrossingPoint(grid, points[0]);
            point1.gridLength = 0.0f;

            int end;
            CellSegment cellSegment = null;
            for(p = 1; p < nPoints; p++)
            {
                point0 = point1;
                point1 = new StsGridCrossingPoint(grid, points[p]);
                point1.computeGridLength(point0);
                ArrayList gridPoints = getGridCrossings(point0, point1);

                if(gridPoints != null)
                {
                    int nGridPoints = gridPoints.size();
                    gridPoint1 = (StsGridCrossingPoint)gridPoints.get(0);
                    for(g = 1; g < nGridPoints; g++)
                    {
                        if(p == 1 && g == 1)
                            end = MINUS;
                        else if(p == nPoints-1 && g == nGridPoints-1)
                            end = PLUS;
                        else
                            end = NONE;

                        gridPoint0 = gridPoint1;
                        gridPoint1 = (StsGridCrossingPoint)gridPoints.get(g);
                        cellSegment = new CellSegment(edge, gridPoint0, gridPoint1, end);
                        addCellSegment(cellSegment);
                    }
                }
            }
            if(point1 != null) intersectedEdge.setGridLength(point1.gridLength);
        }
        catch(Exception e)
        {
            StsException.outputException("StsIntersectSectionEdges.addEdge() failed.",
                e, StsException.WARNING);
        }
    }
/*
    private void addEdge(StsSectionEdge edge)
    {
        StsGridCrossingPoint gridPoint0, gridPoint1;
        float rowF0, colF0, rowF1, colF1;
        int row, col, p, g;

        try
        {
            edge.setPointIndexFs();
            IntersectedEdge intersectedEdge = new IntersectedEdge(edge);
            Integer key = new Integer(edge.index());
            intersectedEdges.put(key, intersectedEdge);

            StsPoint[] points = edge.getPoints();
            int edgeIndex = edge.index();
            int nPoints = points.length;

            ArrayList cellSegments = new ArrayList(100);
            gridPoint1 = new StsGridCrossingPoint(grid, points[0]);

            int end;
            CellSegment cellSegment = null;
            for(p = 1; p < nPoints; p++)
            {
                gridPoint0 = gridPoint1;
                gridPoint1 = new StsGridCrossingPoint(grid, points[p]);
                ArrayList gridPoints = getGridCrossings(gridPoint0, gridPoint1);

                if(gridPoints != null)
                {
                    int nGridPoints = gridPoints.size();
                    gridPoint1 = (StsGridCrossingPoint)gridPoints.get(0);
                    for(g = 1; g < nGridPoints; g++)
                    {
                        if(p == 1 && g == 1)
                            end = MINUS;
                        else if(p == nPoints-1 && g == nGridPoints-1)
                            end = PLUS;
                        else
                            end = NONE;

                        gridPoint0 = gridPoint1;
                        gridPoint1 = (StsGridCrossingPoint)gridPoints.get(g);
                        cellSegment = new CellSegment(edge, gridPoint0, gridPoint1, end);
                        addCellSegment(cellSegment);
                    }
                }
            }

            if(cellSegment != null)
            {
                float gridLength = cellSegment.getTotalGridLength();
                intersectedEdge.setGridLength(gridLength);
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsIntersectSectionEdges.addEdge() failed.",
                e, StsException.WARNING);
        }
    }
*/
    private void addCellSegment(CellSegment cellSegment)
    {
        StsGridCrossingPoint gridPoint0, gridPoint1;
        gridPoint0 = cellSegment.gridPoint0;
        gridPoint1 = cellSegment.gridPoint1;

        int row = (int)Math.min(gridPoint0.iF, gridPoint1.iF);
        row = Math.min(row, nRows-2);
        int col = (int)Math.min(gridPoint0.jF, gridPoint1.jF);
        col = Math.min(col, nCols-2);

//        gridPoint1.gridLength = gridPoint0.gridLength + cellSegment.getDeltaGridLength();

        if(grid.isInsideRowColRange(row, col))
        {
            if(edgeCells[row][col] == null) edgeCells[row][col] = new EdgeCell(row, col);
            edgeCells[row][col].addCellObject(cellSegment);
        }
    }

    private void createIntersectionList()
    {
        intersections = new ArrayList(100);

        for(int row = 0; row < nRows; row++)
            for(int col = 0; col < nCols; col++)
                if(edgeCells[row][col] != null) intersectCellEdges(row, col);
    }

    private void intersectCellEdges(int row, int col)
    {
        Intersection intersection;

        TreeMap extendedIntersections = null;

        EdgeCell edgeCell = edgeCells[row][col];
        Object cellObject;
        CellSegment cellSegment0, cellSegment1;

        Object[] cellObjects = edgeCell.getCellObjects();
        int nObjects = cellObjects.length;
        for(int n0 = 0; n0 < nObjects; n0++)
        {
            cellObject = cellObjects[n0];
            if(!(cellObject instanceof CellSegment)) continue;

            cellSegment0 = (CellSegment)cellObject;
            boolean isExtendable = cellSegment0.isExtendable();

            if(isExtendable) extendedIntersections = new TreeMap(); // a list of extended intersections sorted by gridLength along edge
            for(int n1 = 0; n1 < nObjects; n1++)
            {
                boolean crossOK = n1 > n0;  // don't try to intersect the same pair if CROSSed
                cellObject = cellObjects[n1];
                if(!(cellObject instanceof CellSegment)) continue;
                cellSegment1 = (CellSegment)cellObject;
                if(cellSegment0.edge == cellSegment1.edge) continue;
                intersection = intersectSegments(cellSegment0, cellSegment1);
                if(intersection != null) // if != null, type is CROSS, MINUS, or PLUS
                {
                    if(intersection.type == CROSS)
                    {
                        if(crossOK) intersections.add(intersection);
                    }
                    else if(isExtendable)
                    {
                        Float key = new Float(intersection.gridLengths[0]);
                        extendedIntersections.put(key, intersection);
                    }
                }
            }

            if(isExtendable)
            {
                if(extendedIntersections.size() == 0)
                {
                    EdgeCell[] neighborEdgeCells = getCrossingEdgeCells(cellSegment0, edgeCell);
                    if(neighborEdgeCells == null) continue;

                    int nCells = neighborEdgeCells.length;
                    if(nCells == 0) continue;

                    for(int c = 0; c < nCells; c++)
                    {
                        EdgeCell neighborCell = neighborEdgeCells[c];
                        Object[] neighborCellObjects = neighborCell.getCellObjects();
                        int nNeighborObjects = neighborCellObjects.length;

                        for(int n1 = 0; n1 < nNeighborObjects; n1++)
                        {
                            cellObject = neighborCellObjects[n1];
                            if(!(cellObject instanceof CellSegment)) continue;
                            cellSegment1 = (CellSegment)cellObject;
                            if(cellSegment0.edge == cellSegment1.edge) continue;
                            {
                                intersection = intersectSegments(cellSegment0, cellSegment1);
                                if(intersection != null && intersection.type != CROSS) // type must be MINUS or PLUS
                                {
                                    Float key = new Float(intersection.gridLengths[0]);
                                    extendedIntersections.put(key, intersection);
                                }
                            }
                        }
                    }
                }

                // extendedIntersections are automatically sorted by length, so get
                // either first or last depending on direction

                Object[] extendedObjects = extendedIntersections.values().toArray();
                if(extendedObjects == null) continue;
                int nExtendedObjects = extendedObjects.length;
                if(extendedObjects.length < 1) continue;

                int end = cellSegment0.end;
                if(end == PLUS)
                {
                    intersection = (Intersection)extendedObjects[0];
                    intersections.add(intersection);
                }
                else if(end == MINUS)
                {
                    intersection = (Intersection)extendedObjects[nExtendedObjects-1];
                    intersections.add(intersection);
                }
            }
        }
    }

    EdgeCell[] getCrossingEdgeCells(CellSegment cellSegment, EdgeCell excludeCell)
    {
        float dGridLength = cellSegment.getDeltaGridLength();
        if(dGridLength <= 0.0f) return null;

        float df = maxGridGap/dGridLength;
        float ff;
        ArrayList gridPoints;
        StsGridCrossingPoint gridPointX;
        if(cellSegment.end == MINUS)
        {
            ff = -df;
            gridPointX = new StsGridCrossingPoint(cellSegment.gridPoint0, cellSegment.gridPoint1, ff);
            gridPoints = getGridCrossings(gridPointX, cellSegment.gridPoint0);
        }
        else
        {
            ff = 1.0f + df;
            gridPointX = new StsGridCrossingPoint(cellSegment.gridPoint0, cellSegment.gridPoint1, ff);
            gridPoints = getGridCrossings(cellSegment.gridPoint1, gridPointX);
        }

        if(gridPoints == null) return null;

        int nCells =  gridPoints.size();
        EdgeCell[] crossingCells = new EdgeCell[nCells];
        int lastRow = -1;
        int lastCol = -1;
        int c = 0;
        StsGridCrossingPoint gridPoint;
        for(int n = 0; n < nCells; n++)
        {
            gridPoint = (StsGridCrossingPoint)gridPoints.get(n);
            int[] rowCol = gridPoint.getLowerLeftRowCol();
            int row = rowCol[0];
            int col = rowCol[1];
            if(!insideGrid(row, col))
            {
                // StsException.systemError(this, "getCrossingEdgeCells", "row or col out of range");
                continue;
            }
            EdgeCell edgeCell = edgeCells[row][col];
            if( (row != lastRow || col != lastCol) && edgeCell != null && edgeCell != excludeCell)
                crossingCells[c++] = edgeCells[row][col];
            lastRow = row;
            lastCol = col;
        }
        if(c != nCells)
        {
            EdgeCell[] tempCells = crossingCells;
            crossingCells = new EdgeCell[c];
            System.arraycopy(tempCells, 0, crossingCells, 0, c);
        }
        return crossingCells;
    }

    Intersection intersectSegments(CellSegment segment0, CellSegment segment1)
    {
        float[] rowColF;

        if(segment0.edge == segment1.edge) return null;
        if(edgesAlreadyConnected(segment0, segment1)) return null;

        float[] p00 = new float[] { segment0.gridPoint0.iF, segment0.gridPoint0.jF };
        float[] p01 = new float[] { segment0.gridPoint1.iF, segment0.gridPoint1.jF };
        float[] p10 = new float[] { segment1.gridPoint0.iF, segment1.gridPoint0.jF };
        float[] p11 = new float[] { segment1.gridPoint1.iF, segment1.gridPoint1.jF };

        p01[0] -= p00[0];
        p01[1] -= p00[1];
        p11[0] -= p10[0];
        p11[1] -= p10[1];

        float[] factors = new float[2];

        rowColF = StsMath.linePVIntersectXY(p00, p01, p10, p11, factors);
        if(rowColF == null) return null;

        return intersectionConstructor(rowColF[0], rowColF[1], segment0, segment1, factors);
    }

    private boolean edgesAlreadyConnected(CellSegment segment0, CellSegment segment1)
    {
        int end;
        StsSectionEdge edge0, edge1;
        StsSurfaceVertex vertex;

        edge0 = segment0.edge;
        edge1 = segment1.edge;

        end = segment0.end;
        if(end == MINUS)
        {
            vertex = edge0.getPrevVertex();
            if(vertex.isConnectedToEdge(edge1)) return true;
        }
        else if(end == PLUS)
        {
            vertex = edge0.getNextVertex();
            if(vertex.isConnectedToEdge(edge1)) return true;
        }

        end = segment1.end;
        if(end == MINUS)
        {
            vertex = edge1.getPrevVertex();
            if(vertex.isConnectedToEdge(edge0)) return true;
        }
        else if(end == PLUS)
        {
            vertex = edge1.getNextVertex();
            if(vertex.isConnectedToEdge(edge0)) return true;
        }

        return false;
    }

    private void makeIntersections()
    {
        Iterator iter;
        Intersection intersection;
        IntersectedEdge intersectedEdge;

        int nIntersections = intersections.size();
        if(nIntersections <= 0) return;

        iter = intersections.iterator();
        while(iter.hasNext())
        {
            intersection = (Intersection)iter.next();
            intersection.addToIntersectedEdges();
        }

        iter = intersectedEdges.values().iterator();
        while(iter.hasNext())
        {
            intersectedEdge = (IntersectedEdge)iter.next();
            intersectedEdge.orderIntersections();
        }

        iter = intersectedEdges.values().iterator();
        while(iter.hasNext())
        {
            intersectedEdge = (IntersectedEdge)iter.next();
            intersectedEdge.makeIntersections();
        }

        iter = intersections.iterator();
        while(iter.hasNext())
        {
            intersection = (Intersection)iter.next();
            intersection.connectEdges();
        }
    }

    /** Access method to build a list of gridCrossing points */
    public ArrayList getGridCrossings(StsGridCrossingPoint gridPoint0, StsGridCrossingPoint gridPoint1)
    {
        StsGridCrossings gridCrossing = new StsGridCrossings(gridPoint0, gridPoint1, true);
        return gridCrossing.gridPoints;
    }

    class EdgeCell
    {
        int row;
        int col;
        ArrayList cellObjects = new ArrayList(2);

        EdgeCell(int row, int col)
        {
            this.row = row;
            this.col = col;
        }

        void addCellObject(CellObject cellObject)
        {
            cellObjects.add(cellObject);
        }

        Object[] getCellObjects()
        {
            return cellObjects.toArray();
        }
    }

    class CellObject
    {
        CellObject nextCellObject = null;

        CellObject()
        {
        }

        void addCellObject(CellObject cellObject)
        {
            if(nextCellObject == null)
                nextCellObject = cellObject;
            else
                nextCellObject.addCellObject(cellObject);
        }
    }

    class CellSegment extends CellObject
    {
        StsSectionEdge edge;
        StsGridCrossingPoint gridPoint0;
        StsGridCrossingPoint gridPoint1;
        int end;

        CellSegment(StsSectionEdge edge, StsGridCrossingPoint gridPoint0, StsGridCrossingPoint gridPoint1, int end)
        {
            super();
            this.edge = edge;
            this.gridPoint0 = gridPoint0;
            this.gridPoint1 = gridPoint1;
            this.end = end;
        }

        float getGridLength(float segmentF)
        {
            return gridPoint0.gridLength + segmentF*(gridPoint1.gridLength - gridPoint0.gridLength);
        }

        float getDeltaGridLength()
        {
            float dRow = gridPoint1.iF - gridPoint0.iF;
            float dCol = gridPoint1.jF - gridPoint0.jF;
            return (float)Math.sqrt(dRow*dRow + dCol*dCol);
        }

        float getTotalGridLength()
        {
            return gridPoint1.gridLength;
        }

        StsPoint getPoint(float segmentF)
        {
            StsPoint point = new StsPoint(4);
            point.interpolatePoints(gridPoint0.point, gridPoint1.point, segmentF);
            return point;
        }

        float getGridLengthGap(float factor)
        {
            if(factor < 0.0f) return -factor*getDeltaGridLength();
            else if(factor > 1.0f) return factor*getDeltaGridLength() - 1.0f;
            else return StsParameters.largeFloat;
        }

        boolean isExtendable()
        {
            StsSurfaceVertex vertex;

            if(end == NONE) return false;

            if(end == MINUS)
                vertex = edge.getPrevVertex();
            else
                vertex = edge.getNextVertex();

            return !vertex.isFullyConnected();
        }
    }

    Intersection intersectionConstructor(float rowF, float colF, CellSegment segment0, CellSegment segment1, float[] factors)
    {
        CellSegment[] segments = new CellSegment[] { segment0, segment1 };
        Intersection intersection =  new Intersection(rowF, colF, segments, factors);
        if(intersection.type == NONE) return null;
        else return intersection;
    }

    class Intersection
    {
        float rowF, colF;
        CellSegment[] segments;
        float[] factors;
        float[] gridLengths;
        boolean[] intersected;
        float[][] dEdgeLengths;
        int type = NONE;
        int splitIndex = -1;
        StsSectionEdge[] oldEdges; // these are the original edges crossing at this intersection
        StsSectionEdge[] newEdges; // another intersection may split edges at this intersection: these are the new edges
        StsSectionEdge[] splitEdges; // at this intersection, these are the new edges which are split
        StsSectionEdge unSplitEdge;  // at this intersection, this is the unsplitEdge if CROSS or the intersected edge if not

        Intersection(float rowF, float colF, CellSegment[] segments, float[] factors)
        {
            int end;

            intersected = new boolean[2];
            intersected[0] = StsMath.betweenInclusive(factors[0], 0.0f, 1.0f);
            intersected[1] = StsMath.betweenInclusive(factors[1], 0.0f, 1.0f);

            if(!intersected[0] && !intersected[1]) return; // type == NONE

            this.rowF = rowF;
            this.colF = colF;
            this.segments = segments;
            this.factors = factors;

            gridLengths = new float[2];
            gridLengths[0] = segments[0].getGridLength(factors[0]);
            gridLengths[1] = segments[1].getGridLength(factors[1]);

            oldEdges = new StsSectionEdge[] { segments[0].edge, segments[1].edge };
            newEdges = new StsSectionEdge[] { segments[0].edge, segments[1].edge };

            dEdgeLengths = new float[2][2];

            if(intersected[0] && intersected[1]) type = CROSS;
            else if(intersected[0])
            {
                checkExtendEdge(1);
                unSplitEdge = oldEdges[0];
            }
            else if(intersected[1])
            {
                checkExtendEdge(0);
                unSplitEdge = oldEdges[1];
            }
        }

        private void checkExtendEdge(int extendedEdgeIndex)
        {
            float extendedGridLength;

            int end = segments[extendedEdgeIndex].end;
            float intersectionGridLength = gridLengths[extendedEdgeIndex];
            if(end == PLUS)
            {
                extendedGridLength = intersectionGridLength - segments[extendedEdgeIndex].gridPoint1.gridLength;
                if(extendedGridLength > maxGridGap) return;  // type = NONE;
                type = PLUS;
            }
            else if(end == MINUS)
            {
                extendedGridLength = segments[extendedEdgeIndex].gridPoint0.gridLength - intersectionGridLength;
                if(extendedGridLength > maxGridGap) return;  // type = NONE;
                type = MINUS;
            }
        }

        void addToIntersectedEdges()
        {
            IntersectedEdge intersectedEdge;
            for(int n = 0; n < 2; n++)
            {
                StsSectionEdge edge = oldEdges[n];
                Integer key = new Integer(edge.getIndex());

                intersectedEdge = (IntersectedEdge)intersectedEdges.get(key);
                if(intersectedEdge == null)
                {
                    StsException.systemError("Couldn't find intersectedEdge in intersectedEdges HashMap for edge: " + edge.getLabel());
                    return;
                }
                intersectedEdge.addIntersection(gridLengths[n], this);
            }
        }

        void setEdgeSegmentLengths(StsSectionEdge edge, Object prevLength, Object length, Object nextLength)
        {
            int nEdge = getOldEdgeIndex(edge);
            if(nEdge < 0) return;

            float dPrevLength = ((Float)length).floatValue() - ((Float)prevLength).floatValue();
            float dNextLength = ((Float)nextLength).floatValue() - ((Float)length).floatValue();

            dEdgeLengths[nEdge][0] = dPrevLength;
            dEdgeLengths[nEdge][1] = dNextLength;
        }

        int getOldEdgeIndex(StsSectionEdge edge)
        {
            if(oldEdges[0] == edge) return 0;
            else if(oldEdges[1] == edge) return 1;
            return -1;
        }

        int getOtherEdgeIndex(int index)
        {
            if(index == 0) return 1;
            else if(index == 1) return 0;
            return -1;
        }

        /* Since the intersection appears in intersection lists of both edges, do
           the intersection only if this edge is being split or intersecting the
           other edge.
        */
        StsSectionEdge[] splitEdge(StsSectionEdge edge)
        {
            StsSectionEdge splitEdge;
            if(splitIndex != -1) return null;

            int edgeIndex = getOldEdgeIndex(edge);
            if(edgeIndex < 0) return null;

            int otherEdgeIndex = getOtherEdgeIndex(edgeIndex);

            if(intersected[0] && intersected[1]) // edges cross: find shortest
            {
                float edgeMin = Math.min(dEdgeLengths[edgeIndex][0], dEdgeLengths[edgeIndex][1]);
                float otherEdgeMin = Math.min(dEdgeLengths[otherEdgeIndex][0], dEdgeLengths[otherEdgeIndex][1]);

                if(edgeMin > otherEdgeMin) return null;

                splitEdge = newEdges[edgeIndex];

                StsPoint point = segments[edgeIndex].getPoint(factors[edgeIndex]);
                splitEdges = splitEdge.splitEdge(point);
                if(splitEdges == null) return null;

                unSplitEdge = newEdges[otherEdgeIndex];

                for(int n = 0; n < 2; n++)
                {
                    if(dEdgeLengths[edgeIndex][n] >= maxGridGap) continue;
                    splitEdges[n].delete();
                    splitEdges[n] = null;
                }
                splitIndex = edgeIndex;

                return splitEdges;
            }
            else if(intersected[otherEdgeIndex])
            {
                extendEdge(edgeIndex);
                unSplitEdge = newEdges[otherEdgeIndex]; // this edge is intersected
            }
        /*
            else if(intersected[0])
                extendEdge(1);
            else if(intersected[1])
                extendEdge(0);
        */

            return splitEdges;
        }

        // hack method for now
        StsSectionEdge getIntersectedEdge()
        {
            return unSplitEdge;
        }

        StsSectionEdge getNewEdgeFromOldEdge(StsSectionEdge oldEdge)
        {
            if(oldEdges[0] == oldEdge) return newEdges[0];
            else if(oldEdges[1] == oldEdge) return newEdges[1];
            else
            {
                StsException.systemError("StsIntersectSectionEdges.Intersection.getNewEdgeFromOldEdge() failed." +
                " Couldn't find oldEdge corresponding to: " + oldEdge.getLabel());
                return null;
            }
        }
        // extend edge indicated by extendEdgeIndex against the other edge.
        // If extended in PLUS direction, splitEdges[0] will be this extended edge
        // and if in MINUS direction, splitEdges[1] will be this extend edge; the
        // intersected edge will be the unsplit edge.

        void extendEdge(int extendEdgeIndex)
        {
            try
            {
                int extendEnd = segments[extendEdgeIndex].end;  // edge has to be extended in this direction to intersect
                if(extendEnd == NONE)
                {
                    StsException.systemError("StsIntersectSectionEdges.Intersection.extendEdge() failed. " +
                        "extendEdge == NONE for intersection at rowF: " + rowF + " colF: " + colF);
                    return;
                }

                int intersectEdgeIndex = getOtherEdgeIndex(extendEdgeIndex);
            /*
                float gridLengthGap = segments[extendEdgeIndex].getGridLengthGap(factors[extendEdgeIndex]);
                if(gridLengthGap > maxGridGap) return;

                if(factors[extendEdgeIndex] < 0.0f) extendEnd = MINUS;
                else                                extendEnd = PLUS;
            */
                StsPoint point = segments[intersectEdgeIndex].getPoint(factors[intersectEdgeIndex]);
                StsSectionEdge extendEdge = newEdges[extendEdgeIndex];
                extendEdge.extendEdge(point, extendEnd);
                splitEdges = new StsSectionEdge[2];
                if(extendEnd == PLUS)
                    splitEdges[0] = extendEdge;
                else
                    splitEdges[1] = extendEdge;

                splitIndex = extendEdgeIndex;

                unSplitEdge = newEdges[intersectEdgeIndex];
            }
            catch(Exception e)
            {
                StsException.outputException("StsIntersectSectionEdges.Intersection.extendEdge() failed " +
                    "for intersection at rowF: " + rowF + " colF: " + colF, e, StsException.WARNING);
            }
        }

        // For this intersection there are the two original edges, two new edges
        // assigned when there is a split somewhere else which changes the original
        // edges, and if this intersection is already split, there are two split edges
        // and one unsplit edge.  oldEdge is the original edge.  If the oldEdge is
        // already changed to a newEdge, than set that newEdge to this new newEdge.
        // In addition, if this intersection has been split (splitIndex >= 0), and
        // the oldEdge corresponds to the splitIndex, than newSplitIndex indicates
        // which of the splitEdges needs to be reassigned; otherwise set the unSplitEdge
        // to this newEdge.

        boolean setNewEdge(StsSectionEdge originalEdge, StsSectionEdge oldEdge, StsSectionEdge newEdge)
        {
            int nOldEdge = getOldEdgeIndex(originalEdge);
            if(splitIndex == nOldEdge) // we have reached a split on this edge
            {
                if(splitEdges[0] == oldEdge) splitEdges[0] = newEdge;
                else if(splitEdges[1] == oldEdge) splitEdges[1] = newEdge;
                return true;
            }
            else if(splitIndex != -1) // thus edge splits the crossing edge
                unSplitEdge = newEdge;
            else // this intersection has not been split yet
                newEdges[nOldEdge]= newEdge;

            return false; // we have hit a terminating split
        }

        void connectEdges()
        {
            StsSurfaceVertex vertex;
            if(splitEdges == null) return;
            if(splitEdges[0] != null)
            {
                vertex = splitEdges[0].getNextVertex();
                vertex.setSectionEdge(unSplitEdge);
            }
            if(splitEdges[1] != null)
            {
                vertex = splitEdges[1].getPrevVertex();
                vertex.setSectionEdge(unSplitEdge);
            }
        }
    }

    class IntersectedEdge
    {
        StsSectionEdge edge;
        float gridLength;
        TreeMap intersections; // a list of intersections sorted by gridLength along edge
        Intersection[] intersects;

        IntersectedEdge(StsSectionEdge edge)
        {
            this.edge = edge;
            intersections = new TreeMap();
        }

        public int hashCode()
        {
            return edge.getIndex();
        }

        public boolean equals(Object other)
        {
            if(other.getClass() != IntersectedEdge.class) return false;
            return hashCode() == ((IntersectedEdge)other).hashCode();
        }

        void addIntersection(float gridLength, Intersection intersection)
        {
            Float key = new Float(gridLength);
            intersections.put(key, intersection);
        }

        void setGridLength(float gridLength)
        {
            this.gridLength = gridLength;
        }

        void orderIntersections()
        {
            Intersection intersection;

            int nIntersections = intersections.size();
            if(nIntersections == 0) return;

            intersects = new Intersection[nIntersections];

            Iterator iter = intersections.values().iterator();
            int n = 0;
            while(iter.hasNext())
                intersects[n++] = (Intersection)iter.next();

            Object[] lengths = intersections.keySet().toArray();

            Float startLength = new Float(0.0f);
            Float endLength = new Float(gridLength);

            for(n = 0; n < nIntersections; n++)
                intersects[n].setEdgeSegmentLengths(edge, startLength, lengths[n], endLength);
        /*
            Object prevLength;
            Object length = (Object)startLength;
            Object nextLength = lengths[0];
            for(n = 0; n < nIntersections; n++)
            {
                prevLength = length;
                length = nextLength;

                if(n < nIntersections-1)
                    nextLength = lengths[n+1];
                else
                    nextLength = endLength;

                intersection = intersects[n];
                intersection.setEdgeSegmentLengths(edge, prevLength, length, nextLength);
            }
        */
        }

        void makeIntersections()
        {
            if(intersects == null) return;

            int nIntersections = intersects.length;
            if(nIntersections > 2)
            {
                StsException.systemError(this, "makeIntersections", "This sectionEdge " + edge.toString() +
                    " has more than two intersections. Construct not implemented yet.");
                return;
            }
            for(Intersection intersect : intersects)
            {
                if(intersect.type == CROSS)
                {
                    StsException.systemError(this, "makeIntersections", "This sectionEdge " + edge.toString() +
                        " has a crossing intersection. Construct not implemented yet.");
                    return;
                }
            }
            // Assumption at this point is that we have one or two intersections which involve simply extending the
            // existing section in the minus and/or the plus direction(s).
            // So for now, we will ignore the splitEdges and simply extend the ends of the existing section as needed.
            StsSection section = edge.getSection();
            if(section == null)
            {
                    StsException.systemError(this, "makeIntersections", "This sectionEdge " + edge.toString() +
                        " has no section. Construct not implemented yet.");
                    return;
            }
            boolean sectionExtended = false;
            // we currently have the same intersect object on both intersecting edges
            // for a "T" intersection, we ignore the intersection on the "bar" edge and attach the "vertical" edge to the "bar" edge
            for(Intersection intersect : intersects)
            {
                StsSectionEdge intersectedEdge = intersect.getIntersectedEdge();
                if(edge == intersectedEdge) continue; // ignore the intersection if the edge && intersected edge are the same
                StsSection intersectedSection = intersectedEdge.getSection();
                if(intersectedSection == null)
                {
                     StsException.systemError(this, "makeIntersections", "This sectionEdge " + edge.toString() +
                        " has no intersected section.");
                    continue;
                }
                StsLine projectedLine;
                if(intersect.type == MINUS)
                    projectedLine = (StsLine)section.getLines().getFirst();
                else
                    projectedLine = (StsLine)section.getLines().getLast();
                projectedLine.setOnSection(intersectedSection);
                projectedLine.projectVerticesToSection();
                // projectedLine.projectToSection();
                section.addLineToSectionSide(intersect.type);
                sectionExtended = true;
            }
            if(sectionExtended)
            {
                section.reinitialize();
                //model.addMethodCmd(section, "reinitialize");
            }
            /*
            for(int n = 0; n < nIntersections; n++)
            {
                StsSectionEdge[] splitEdges = intersects[n].splitEdge(edge);
                if(splitEdges != null)
                {
                    StsSectionEdge preSplitEdge = intersects[n].getNewEdgeFromOldEdge(edge);
                    if(preSplitEdge == null) continue;

                    if(splitEdges[0] != null)
                    {
                        for(int nn = n-1; nn >= 0; nn--)
                            if(intersects[nn].setNewEdge(edge, preSplitEdge, splitEdges[0])) break;
                    }

                    if(splitEdges[1] != null)
                    {
                        for(int nn = n+1; nn < nIntersections; nn++)
                            if(intersects[nn].setNewEdge(edge, preSplitEdge, splitEdges[1])) break;
                    }
                }
            }
            */
        }
    }
}
