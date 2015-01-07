package com.Sts.DBTypes;

import com.Sts.DB.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.View3d.*;
import com.Sts.SeismicAttributes.*;
import com.Sts.Types.*;
import com.Sts.Utilities.Seismic.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jun 12, 2009
 * Time: 11:17:20 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsPatchGrid extends StsXYGridBoundingBox implements Comparable<StsPatchGrid>, StsSerializable, StsXYSurfaceLinkGridable
{
    private static final long serialVersionUID = 1L;
    /** this is the parent id of this child patchGrid */
    int id;
    // /** the parent multilayered patch grid is divided into single-valued sheets; each sheet is divided into child patchGrids; this is the sheet number for this child grid */
    // int nSheet;
    /** the final id of this patch in gridList; this patch may be a layer from an original patch whose id is id */
    int idFinal;
    /** type of this grid: Min, Max +Zero-crossing, or -Zero-crossing */
    byte patchType;
    //TODO make patchPoints a row&col sorted list for compactness; add methods to compare patchPoint lists for occupancy, etc.
    StsPatchPoint[][] patchPoints;
    boolean isVisible = true;
    StsPatchVolume patchVolume;
    int nPatchPoints;
    int nRows;
    int nCols;
    float dataMin;
    float dataMax;
    float zMin = StsParameters.largeFloat;
    float zMax = -StsParameters.largeFloat;
    /** included for debugging purposes only */
    transient int nSheets;
    transient float[][] curvature;
    transient float[][] pointsZ;

    transient StsDiamondStrips diamondStrips;

    transient public int nValuePatchPoints = 0;
    transient public double sum = 0;
    /** points on this grid are in a hashMap with a hash key whose value is the grid index: col + row*patchVolume.nCols */
    // transient HashMap<Integer, StsPatchPoint> patchPointsHashMap;
    transient ArrayList<StsPatchConnection> connectionsList;
    transient HashMap<Integer, StsPatchPoint> pointHashMap;
    transient StsList tStrips;
    transient float[][][] tStripNormals;
    /** max size of all patches generated; printed out for diagnostics */
    static int maxGridSize = 0;
    /** index to be assigned to next patchGrid created; this is the index during construction for a multilayered patch;
     *  incremented for each one; reset on initialization
     */
    static int nextPatchID = 0;
    /** final index in gridList reset on initialization */
    static int nextFinalPatchID = 0;

    static final float nullValue = StsParameters.nullValue;

    // Multiply # pts in SVD to get ChiSqr limit
    static final private double chiSqrMultiplyer = 2;
    //static final private double stdDevFactor = 1.5;
    static final byte FILTER_NONE = 0;
    static final byte FILTER_ON_CHI_SQ = 1;
    static final byte FILTER_ON_STD_DEV = 2;

    static public byte filterType = FILTER_ON_CHI_SQ;

    static public final float badCurvature =  StsQuadraticCurvature.badCurvature;
    static public final float curvatureTest = StsQuadraticCurvature.curvatureTest;

    static boolean sortRowFirst = true;

    static int nOverlappedPoints = 0;
    static int nVolumeRows;
    static int nVolumeCols;
    static ArrayList<OverlapPoint> overlapPoints;

    static final boolean debug = false;
    /** various debug print of patches in rowGrid and prevRowGrid arrays; check if this is not -1 and prints if id matches this value.  Make this -1 if you want no debug prints. */
    static final int debugPatchInitialID = -1;
    /** debugPatchID may change if orginal patch is merged into a new one; when this occurs, set debugCurrentPatchID to new one and track it */
    static int debugPatchID = debugPatchInitialID;

    static final int debugLayerPointRow = -1;
    static final int debugLayerPointCol = -1;

    public StsPatchGrid()
    {
    }
    
    private StsPatchGrid(StsPatchVolume patchVolume, byte patchType)
    {
        this.patchVolume = patchVolume;
        this.patchType = patchType;
    }

    static public void staticInitialize(int nVolumeRows_, int nVolumeCols_)
    {
        maxGridSize = 0;
        nextPatchID = 0;
        nextFinalPatchID = 0;
        nOverlappedPoints = 0;
        nVolumeRows = nVolumeRows_;
        nVolumeCols = nVolumeCols_;
        debugPatchID = debugPatchInitialID;
        if(debug) overlapPoints = new ArrayList();
    }

    static public StsPatchGrid construct(StsPatchVolume patchVolume, byte patchType)
    {
        StsPatchGrid patchGrid = new StsPatchGrid(patchVolume, patchType);
        patchGrid.setup();
        return patchGrid;
    }

    static public StsPatchGrid constructFinal(StsPatchVolume patchVolume, byte patchType, int id, int rowMin, int rowMax, int colMin, int colMax)
    {
        StsPatchGrid patchGrid = new StsPatchGrid(patchVolume, patchType);
        patchGrid.setupFinal(id, rowMin, rowMax, colMin, colMax);
        return patchGrid;
    }

    private void setup()
    {
        id = nextPatchID++;
        // patchPointsHashMap = new HashMap<Integer, StsPatchPoint>();
        connectionsList = new ArrayList<StsPatchConnection>();
        if(debugPatchID != -1 && id == debugPatchID)
            StsException.systemDebug(this, "setup", "patch " + id + " initialized");
    }
    private void setupFinal(int id, int rowMin, int rowMax, int colMin, int colMax)
    {
        this.id = id;
        this.rowMin = rowMin;
        this.rowMax = rowMax;
        this.colMin = colMin;
        this.colMax = colMax;
        nRows = rowMax - rowMin + 1;
        nCols = colMax - colMin + 1;
        patchPoints = new StsPatchPoint[nRows][nCols];
        idFinal = nextFinalPatchID++;
        if(debugPatchID != -1 && id == debugPatchID)
            StsException.systemDebug(this, "setup", "patch " + id + " initialized");
    }

    /** add this connection to the initial (possibly multilayered) patch */
    void addInitialPatchConnection(StsPatchConnection connection)
    {
        addPatchConnection(connection, id);
    }

   /** Add this connection between newPoint and otherPoint to the patchGrid.
     *  Add the two connected points to this patch unless they already belong to another patch.
     *  Add the connection to this patch.
     * @param connection  connection between two points
     */
   void addPatchConnection(StsPatchConnection connection, int id)
    {
        StsPatchPoint newPoint = connection.newPoint;
        StsPatchPoint otherPoint = connection.otherPoint;
        newPoint.patchID = id;
        otherPoint.patchID = id;

        // otherPoint is always before new point in either row or column
        if(otherPoint.row < rowMin) rowMin = otherPoint.row;
        if(otherPoint.col < colMin) colMin = otherPoint.col;
        if(newPoint.row > rowMax) rowMax = newPoint.row;
        if(newPoint.col > colMax) colMax = newPoint.col;

        if(debugPatchID != -1 && (connection.newPoint.patchID == debugPatchID || connection.newPoint.patchID == debugPatchID))
            StsException.systemDebug(this, "addPatchConnection", "connection: " + connection.toString());
        connectionsList.add(connection);

    }

    /**
     * Merge another patchGrid into this by creating a union of the two bounding boxes and
     * copy the floats from the originalGrid and otherGrid into the new one.
     * reset the patchID for each of the merged patchPoints
     *
     * @param mergedPatchGrid the other patch grid to be add to this one.
     */

    public void mergePatchGrid(StsPatchGrid mergedPatchGrid, StsPatchConnection newConnection)
    {
        if(debugPatchID != -1 && (id == debugPatchID || mergedPatchGrid.id == debugPatchID))
        {
            if(mergedPatchGrid.id == debugPatchID)
            {
                StsException.systemDebug(this, "mergePatchGrid", "merged patch " + mergedPatchGrid.toString() + " into patch" + toString() + " will remove patch: " + mergedPatchGrid.id + "\n" +
                    "     Switching debugPatchID to new patch: " + id);
                debugPatchID = id;
            }
            else
                StsException.systemDebug(this, "mergePatchGrid", "merged patch " + mergedPatchGrid.toString() + " into patch" + toString() + " will remove patch: " + mergedPatchGrid.id);
        }
 
        for (StsPatchConnection connection : mergedPatchGrid.connectionsList)
            addInitialPatchConnection(connection);
        addInitialPatchConnection(newConnection);
        mergedPatchGrid.clear();
    }

    public void clear()
    {
        patchPoints = null;
        if(debugPatchID != -1 && id == debugPatchID)
            StsException.systemDebug(this, "clear", "clearing patch: " + toString());
    }

    boolean isDisconnected(int row)
    {
        return rowMax < row;
    }

    boolean isOnePoint()
    {
        return rowMax - rowMin <= 0 && colMax - colMin <= 0;
    }

    boolean isTooSmall(int minNPoints)
    {
        return nPatchPoints < minNPoints;
    }

    public void initializeGrid()
    {
        StsPatchPoint point;

        nRows = rowMax - rowMin + 1;
        nCols = colMax - colMin + 1;
        boolean[][] hasPoint = new boolean[nRows][nCols];
        for (StsPatchConnection connection : connectionsList)
        {
            point = connection.otherPoint;
            hasPoint[point.row - rowMin][point.col - colMin] = true;
            point = connection.newPoint;
            hasPoint[point.row - rowMin][point.col - colMin] = true;
        }
        nPatchPoints = 0;
        for(int row = 0; row < nRows; row++)
            for(int col = 0; col < nCols; col++)
                if(hasPoint[row][col]) nPatchPoints++;
    }

    public void resetIndex(int index)
    {
        if(debugPatchID != -1 && id == debugPatchID)
            StsException.systemDebug(this, "resetIndex", "debugPatch final ID being reset from " + idFinal + " to " + index);
        idFinal = index;
        for(int row = 0; row < nRows; row++)
            for(int col = 0; col < nCols; col++)
                if(patchPoints[row][col] != null)
                    patchPoints[row][col].patchID = idFinal;
    }

    /* This patch contains a connection list where a connection is between a newPoint and the otherPoint.
     * Sort all the connections into a series of sheets where the first sheet contains all first connections,
     *  second sheet has second connections, etc. Each sheet is sorted into row-col order.
     *  All connections within a sheet may not be actually interconnected, so now extract a set of layers
     *  from each sheet where all connections within a layer are interconnected.  Some connections in a sheet
     *  may in fact belong to a layer initially constructed from a previous sheet.
     *  Make successive passes thru the connections assigning all connections which connect to each other
     *  to a layer.  Add this new layer to the gridList.
     *  For each layer, compute the patchPoints array from the connections and delete the connections.
     * @param gridList
     */
    void constructGridArray(ArrayList<StsPatchGrid> gridList)
    {
        if(debugPatchID != -1 && id == debugPatchID)
            StsException.systemDebug(this, "constructGridArray", "patch: " + id);
        patchVolume.nParentGrids++;

        maxGridSize = StsMath.max3(maxGridSize, nRows, nCols);

        try
        {
            // first make a sorted list by row, column, and z
            Collections.sort(connectionsList);
            // assign a sheet number for debugging purposes
            setTraceOrder();
            // clear the patchID numbers from all the patchPoints
            clearPatchPointIDs();
            // local array of patchGrids built from this multilayered patchGrid
            ArrayList<StsPatchGrid> layerGrids = new ArrayList<>();
            int parentID = id;
            StsPatchConnection[] connections = connectionsList.toArray(new StsPatchConnection[0]);
            StsPatchConnection[] pointConnections;
            int nConnections = connections.length;
            StsPatchConnection connection, nextConnection = null;
            // get the connection and the nextConnection (the one after that)
            // if both from same point, increment counter and add both connections using patchID info
            // otherwise just add the one connection
            for(int n = 0; n < nConnections; n++)
            {
                StsPatchGrid layerGrid;

                connection = connections[n];

                if(n+1 < nConnections)
                    nextConnection = connections[n+1];
                else
                    nextConnection = null;

                if(nextConnection != null)
                {
                    if (connection.newPoint != nextConnection.newPoint)
                        nextConnection = null;
                }
                if(nextConnection == null)
                {
                    int id = connection.otherPoint.patchID;
                    StsPatchPoint newPoint = connection.newPoint;
                    if(id == -1) // connection has no associated layer; create one and add both points to it
                    {
                        layerGrid = StsPatchGrid.constructFinal(patchVolume, patchType, parentID, rowMin, rowMax, colMin, colMax);
                        layerGrids.add(layerGrid);
                        layerGrid.addPoint(connection.otherPoint, true);
                        layerGrid.addPoint(newPoint, true);
                        connection.addCorrelation();
                    }
                    else // connection belongs to existing layer; add newPoint to it
                    {
                        layerGrid = getLayerGrid(layerGrids, id);
                        layerGrid.addPoint(layerGrid, layerGrids, newPoint, parentID);
                        connection.addCorrelation();
                    }
                }
                else // nextConnection != null :  we have two connections to the same newPoint
                {
                    n++;
                    StsPatchPoint newPoint = connection.newPoint;  // connection2 is and must be connected to the same point
                    int id1 = connection.otherPoint.patchID;
                    int id2 = nextConnection.otherPoint.patchID;
                    StsPatchGrid layerGrid1, layerGrid2;

                    if(id1 == id2) // both connected points belong to same grid or both to no grid
                    {
                        if(id1 == -1) // need to add all 3 points to a new layer
                        {
                            layerGrid2 = StsPatchGrid.constructFinal(patchVolume, patchType, parentID, rowMin, rowMax, colMin, colMax);
                            layerGrids.add(layerGrid2);
                            addPoint(layerGrid2, layerGrids, connection.otherPoint, parentID);
                            connection.addCorrelation();
                            addPoint(layerGrid2, layerGrids, nextConnection.otherPoint, parentID);
                            nextConnection.addCorrelation();
                            addPoint(layerGrid2, layerGrids, newPoint, parentID);
                        }
                        else // both otherPoints belong to the same grid: add newPoint to it
                        {
                            layerGrid = getLayerGrid(layerGrids, id1);
                            addPoint(layerGrid, layerGrids, newPoint, parentID);
                            connection.addCorrelation();
                            nextConnection.addCorrelation();
                        }
                    }
                    else if(id1 == -1)  // id2 != -1 : nextConnection.otherPoint belongs to a layerGrid2; add the other 2 points to it
                    {
                        layerGrid2 = getLayerGrid(layerGrids, id2);
                        addPoint(layerGrid2, layerGrids, connection.otherPoint, parentID);
                        connection.addCorrelation();
                        nextConnection.addCorrelation();
                        addPoint(layerGrid2, layerGrids, newPoint, parentID);
                    }
                    else if(id2 == -1) // connection.otherPoint belongs to a layer; add the other 2 points to it
                    {
                        layerGrid1 = getLayerGrid(layerGrids, id1);
                        addPoint(layerGrid1, layerGrids, nextConnection.otherPoint, parentID);
                        connection.addCorrelation();
                        nextConnection.addCorrelation();
                        addPoint(layerGrid1, layerGrids, newPoint, parentID);
                    }
                    else // connections are to two different layerGrids; merge if they don't overlap; otherwise clone point and add to the second grid
                    {
                        layerGrid1 = getLayerGrid(layerGrids, id1);
                        layerGrid2 = getLayerGrid(layerGrids, id2);
                        // If newPoint overlaps point in one of the two grids, add it to the other and don't merge grids; clone connected point from overlapped grid
                        // and add it to the non-overlapped grid.
                        // If newPoint doesn't overlap either grid and the two grids don't overlap, we can merge the two grids; add the newPoint to one and merge them.
                        // If newPoint overlaps both grids, then we need to make a new grid containing just this point with clones of the two connected otherPoints.

                        // overlap flag indicates whether there the newPoint overlaps zero grids (0), layerGrid1 only (1), layerGrid2 only (2), or both (3)
                        int overlapFlag = pointOverlapsGrid(newPoint, layerGrid1, layerGrid2);
                        if(overlapFlag == 0) // newPoint doesn't overlap either grid: try to merge
                        {
                            if(mergeGrids(layerGrid1, layerGrid2, layerGrids))
                            {
                                connection.addCorrelation();
                                addPoint(layerGrid1, layerGrids, newPoint, parentID);
                                nextConnection.addCorrelation();
                            }
                            else // can't merge, arbitrarily add newPoint to layerGrid1; clone connected point in layerGrid2 and add to layerGrid1 also
                            {
                                connection.addCorrelation();
                                addPoint(layerGrid1, layerGrids, newPoint, parentID);
                                nextConnection.otherPoint = nextConnection.otherPoint.cloneClearCorrels();
                                nextConnection.addCorrelation();
                                addPoint(layerGrid1, layerGrids, nextConnection.otherPoint, parentID);
                            }
                        }
                        else if(overlapFlag == 1) // newPoint overlaps layerGrid1, so add it to layerGrid2 and clone layerGrid1 connected point and add to layerGrid2
                        {
                            nextConnection.addCorrelation();
                            addPoint(layerGrid2, layerGrids, newPoint, parentID);
                            connection.otherPoint = connection.otherPoint.cloneClearCorrels();
                            connection.addCorrelation();
                            addPoint(layerGrid2, layerGrids, connection.otherPoint, parentID);
                        }
                        else if(overlapFlag == 2) // newPoint overlaps layerGrid2, so add it to layerGrid1 and clone layerGrid2 connected point and add to layerGrid1
                        {
                            connection.addCorrelation();
                            addPoint(layerGrid1, layerGrids, newPoint, parentID);
                            nextConnection.otherPoint = nextConnection.otherPoint.cloneClearCorrels();
                            nextConnection.addCorrelation();
                            addPoint(layerGrid1, layerGrids, nextConnection.otherPoint, parentID);
                        }
                        else // overlapFlag == 3 : newPoint overlaps both grids, so create new layerGrid, and newPoint and cloned otherPoints to new layerGrid
                        {
                            layerGrid = StsPatchGrid.constructFinal(patchVolume, patchType, parentID, rowMin, rowMax, colMin, colMax);
                            layerGrids.add(layerGrid);
                            addPoint(layerGrid, layerGrids, newPoint, parentID);
                            connection.otherPoint = connection.otherPoint.cloneClearCorrels();
                            connection.addCorrelation();
                            addPoint(layerGrid, layerGrids, connection.otherPoint, parentID);
                            nextConnection.otherPoint = nextConnection.otherPoint.cloneClearCorrels();
                            nextConnection.addCorrelation();
                            addPoint(layerGrid, layerGrids, nextConnection.otherPoint, parentID);
                        }
                    }
                }
            }
            connectionsList = null;

            for(StsPatchGrid grid : layerGrids)
            {
                if(grid.nPatchPoints == 0) return;
                grid.resize();
                gridList.add(grid);
                if(debugPatchID != -1 && id == debugPatchID)
                    StsException.systemDebug(this, "constructGridArray", "patchID" + id + " layer patch id " + grid.idFinal);
            }
            if(debug) StsException.systemDebug(this, "constructGridArray", "Constructed " + layerGrids.size() + " new patch grids from " + nSheets + " sheets.");
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "constructGridArray", toString(), e);
        }
    }

    private int pointOverlapsGrid(StsPatchPoint newPoint, StsPatchGrid layerGrid1, StsPatchGrid layerGrid2)
    {
        int row = newPoint.row - rowMin;
        int col = newPoint.col - colMin;
        StsPatchPoint patchPoint1 = layerGrid1.patchPoints[row][col];
        StsPatchPoint patchPoint2 = layerGrid2.patchPoints[row][col];
        if(patchPoint1 == null)
        {
            if(patchPoint2 == null)
                return 0;
            else
                return 2;
        }
        else if(patchPoint2 == null) // patchPoint1 != null
        {
            return 1;
        }
        else // patchPoint1 != null && patchPoint2 != nul
            return 3;
    }

    /** if both grids have points at the same location which aren't the same, then we can't merge.
     *  If we can merge, then add the points from the second grid to the first, add the new connection, and remove the second grid
     * @param layerGrid1 first new child grid
     * @param layerGrid2 second new child grid
     * @param layerGrids the set of child grids being built from the parent grid
     * @return
     */
    private boolean mergeGrids(StsPatchGrid layerGrid1, StsPatchGrid layerGrid2, ArrayList<StsPatchGrid> layerGrids)
    {
        for(int row = 0; row < nRows; row++)
            for(int col = 0; col < nCols; col++)
            {
                StsPatchPoint patchPoint1 = layerGrid1.patchPoints[row][col];
                StsPatchPoint patchPoint2 = layerGrid2.patchPoints[row][col];
                if(patchPoint1 != null && patchPoint2 != null && patchPoint1 != patchPoint2) return false;
            }
       for(int row = 0; row < nRows; row++)
            for(int col = 0; col < nCols; col++)
            {
                StsPatchPoint patchPoint1 = layerGrid1.patchPoints[row][col];
                StsPatchPoint patchPoint2 = layerGrid2.patchPoints[row][col];
                if(patchPoint2 != null)
                {
                    patchPoint2.patchID = layerGrid1.idFinal;
                    if(patchPoint1 != null)
                        mergePatchPoints(patchPoint1, patchPoint2);
                    else
                        layerGrid1.patchPoints[row][col] = patchPoint2;
                }
            }
        layerGrids.remove(layerGrid2);
        return true;
    }

    private void mergePatchPoints(StsPatchPoint patchPoint1, StsPatchPoint patchPoint2)
    {
        if(patchPoint2.rowCorrel != 0.0f)
            patchPoint1.rowCorrel = patchPoint2.rowCorrel;
        if(patchPoint2.colCorrel != 0.0f)
            patchPoint1.colCorrel = patchPoint2.colCorrel;
    }

    private boolean addPoint(StsPatchGrid layerGrid, ArrayList<StsPatchGrid> layerGrids, StsPatchPoint point, int parentID)
    {
        if(layerGrid.addPoint(point, false))
            return true;
        layerGrid = StsPatchGrid.constructFinal(patchVolume, patchType, parentID, rowMin, rowMax, colMin, colMax);
        layerGrids.add(layerGrid);
        return layerGrid.addPoint(point, false);
    }

    private boolean addPoint(StsPatchPoint point, boolean check)
    {
        nPatchPoints++;
        StsPatchPoint currentPoint = patchPoints[point.row-rowMin][point.col-colMin];
        if(currentPoint != null)
        {
            if(check)
                StsException.systemError(this, "addPoint", "Failed. Point already exists:  " + currentPoint.toString());
            return false;
        }
        patchPoints[point.row-rowMin][point.col-colMin] = point;
        point.patchID = idFinal;
        return true;
    }

    private void resize()
    {
        int newRowMin, newRowMax, newColMin, newColMax;
        for(newRowMin = rowMin; newRowMin <= rowMax; newRowMin++)
            if(rowHasNonNull(newRowMin)) break;
        for(newRowMax = rowMax; newRowMax >= rowMin; newRowMax--)
            if(rowHasNonNull(newRowMax)) break;
        for(newColMin = colMin; newColMin <= colMax; newColMin++)
            if(colHasNonNull(newColMin)) break;
        for(newColMax = colMax; newColMax >= colMin; newColMax--)
            if(colHasNonNull(newColMax)) break;
        nRows = newRowMax - newRowMin + 1;
        nCols = newColMax - newColMin + 1;
        StsPatchPoint[][] newPatchPoints = new StsPatchPoint[nRows][nCols];
        try
        {
            for(int row = newRowMin; row <= newRowMax; row++)
            {
                StsPatchPoint[] rowPoints = new StsPatchPoint[nCols];
                System.arraycopy(patchPoints[row-rowMin], newColMin-colMin, rowPoints, 0, nCols);
                newPatchPoints[row-newRowMin] = rowPoints;
            }
        }
        catch(Exception e)
        {
            StsException.systemError(this, "resize");
            return;
        }
        rowMin = newRowMin;
        rowMax = newRowMax;
        colMin = newColMin;
        colMax = newColMax;
        patchPoints = newPatchPoints;
    }

    private boolean rowHasNonNull(int row)
    {
        for(int col = colMin; col <= colMax; col++)
            if(patchPoints[row-rowMin][col-colMin] != null) return true;
        return false;
    }

    private boolean colHasNonNull(int col)
    {
        for(int row = rowMin; row <= rowMax; row++)
            if(patchPoints[row-rowMin][col-colMin] != null) return true;
        return false;
    }

    private void addOtherPoint(StsPatchConnection connection)
    {
        addPoint(connection.otherPoint, true);
        connection.addCorrelation();
    }

    private StsPatchGrid getLayerGrid(ArrayList<StsPatchGrid> layerGrids, int id)
        {
            for(StsPatchGrid layerGrid : layerGrids)
                if(layerGrid.idFinal == id) return layerGrid;
            StsException.systemError(this, "getLayerGrid", "Failed to find layerGrid for id " + id);
            return null;
        }

    private void setTraceOrder()
    {
        int nSheet = 0;
        nSheets = 0;
        StsPatchConnection nextConnection = null;
        for(StsPatchConnection connection : connectionsList)
        {
            StsPatchConnection prevConnection = nextConnection;
            nextConnection = connection;
            if(!nextConnection.sameRowCol(prevConnection))
                nSheet = 0;
            connection.nSheet = nSheet++;
            nSheets = Math.max(nSheets, nSheet);
        }
    }

    private void clearPatchPointIDs()
    {
        for(StsPatchConnection connection : connectionsList)
            connection.clearPointPatchIDs();
    }

    /** get nearest patch whose z at this x,y is just above the z slice plane */
    public float getZDistance(int volumeRow, int volumeCol, float z)
    {
        int volumeRowMin = Math.max(rowMin, volumeRow - 5);
        int volumeRowMax = Math.min(nRows-1 + rowMin, volumeRow + 5);
        int volumeColMin = Math.max(colMin, volumeCol - 5);
        int volumeColMax = Math.min(nCols-1 + colMin, volumeCol + 5);
        float dz = StsParameters.largeFloat;
        for(volumeRow = volumeRowMin; volumeRow <= volumeRowMax; volumeRow++)
        {
            for(volumeCol = volumeColMin; volumeCol <= volumeColMax; volumeCol++)
            {
                float zPatch = getPointZ(volumeRow, volumeCol);
                if(zPatch == nullValue) continue;
                float dzPatch = z - zPatch;
                if(dzPatch < 0.0f) continue;
                if(dzPatch < dz) dz = dzPatch;
            }
        }
        return dz;
    }

    static public void printOverlapPoints()
    {
        for(OverlapPoint overlapPoint : overlapPoints)
            overlapPoint.print();
    }
    class OverlapPoint
    {
        int row, col, id;

        OverlapPoint(int row, int col, int id)
        {
            this.row = row;
            this.col = col;
            this.id = id;
        }

        void print()
        {
            System.out.println("patchPoint overlap at row: " + row + " col: " + col + " id: " + id);
        }
    }

    /** sort first by rowMin and then by colMin. Return 1 if this rowMin&colMin come after other; 0 if equal; -1 otherwise */
    public int compareTo(StsPatchGrid otherGrid)
    {
        if (sortRowFirst)
        {
            if (rowMin > otherGrid.rowMin) return 1;
            if (rowMin < otherGrid.rowMin) return -1;
            // on the same row
            if (colMin > otherGrid.colMin) return 1;
            if (colMin < otherGrid.colMin) return -1;
            return 0;
        }
        else
        {
            if (colMin > otherGrid.colMin) return 1;
            if (colMin < otherGrid.colMin) return -1;
            // on the same col
            if (rowMin > otherGrid.rowMin) return 1;
            if (rowMin < otherGrid.rowMin) return -1;
            return 0;
        }
    }

    public float[][] getPointsZ() 
    {
        if(pointsZ != null) return pointsZ;
        pointsZ = new float[nRows][nCols];
        for(int row = 0; row < nRows; row++)
        {
            for(int col = 0; col < nCols; col++)
            {               
                StsPatchPoint patchPoint = patchPoints[row][col];
                if(patchPoint == null)
                    pointsZ[row][col] = nullValue;
                else
                    pointsZ[row][col] = patchPoints[row][col].z;
            }
        }
        return pointsZ;
    }

    public int getGridSize() { return nRows*nCols; }

    public int[] getGridPointsUsed()
    {
        int nUsed = 0;
        int nActualUsed = 0;
        for(int row = 0; row < nRows; row++)
        {
            int colStart = 0, colEnd = nCols - 1;
            for (int col = 0; col < nCols; col++)
            {
                if (patchPoints[row][col] != null)
                {
                    colStart = col;
                    break;
                }
            }
            for (int col = nCols - 1; col > 0; col--)
            {
                if (patchPoints[row][col] != null)
                {
                    colEnd = col;
                    break;
                }
            }

            for(int col = colStart; col <= colEnd; col++)
                if(patchPoints[row][col] != null) nActualUsed++;

            nUsed += colEnd - colStart + 1;
        }
        return new int[] { nUsed, nActualUsed };
    }

    public float getDataMin()
    {
        return dataMin;
    }
    public float getDataMax()
    {
        return dataMax;
    }
    public float getzMax()
    {
        return zMax;
    }
    public float getzMin()
    {
        return zMin;
    }

    public boolean computeCurvature(float xInc, float yInc, byte curveType, int filterSize, int minNPoints)
    {
        dataMin = StsPatchVolume.largeFloat;
        dataMax = -StsPatchVolume.largeFloat;

        curvature = new float[nRows][nCols];
        for(int row = 0; row < nRows; row++)
            Arrays.fill(curvature[row], nullValue);

        if(nPatchPoints < minNPoints) return false;

        int halfWindow = filterSize/2;
        // Determine quadratic coefficients for this neighborhood

        float[][] fitPoints = new float[filterSize*filterSize][3];
        for (int volumeRow = rowMin; volumeRow <= rowMax; volumeRow++)
        {
            for (int volumeCol = colMin; volumeCol <= colMax; volumeCol++)
            {
                int nFitPoints = 0;  // number of equations
                int patchPointRow = volumeRow - rowMin;
                int patchPointCol = volumeCol - colMin;
                float zc = getPointZ(volumeRow, volumeCol);
                if (zc == StsParameters.nullValue) continue;
                int patchRowMin = Math.max(0, patchPointRow - halfWindow);
                int patchRowMax = Math.min(nRows-1, patchPointRow + halfWindow);
                int patchColMin = Math.max(0, patchPointCol - halfWindow);
                int patchColMax = Math.min(nCols-1, patchPointCol + halfWindow);
                float y = (patchRowMin - patchPointRow)*yInc;
                for (int patchRow = patchRowMin; patchRow <= patchRowMax; patchRow++, y += yInc)
                {
                    float x = (patchColMin - patchPointCol)*xInc;
                    for (int patchCol = patchColMin; patchCol <= patchColMax; patchCol++, x += xInc)
                    {
                        StsPatchPoint patchPoint = patchPoints[patchRow][patchCol];
                        if(patchPoint == null) continue;
                        float z = patchPoint.z;
                        if (z == StsParameters.nullValue) continue;
                        fitPoints[nFitPoints][0] = x;
                        fitPoints[nFitPoints][1] = y;
                        fitPoints[nFitPoints][2] = z - zc;
                        nFitPoints++;
                    }
                }
                if (nFitPoints < minNPoints) continue;

                if(!StsQuadraticCurvature.computeSVD(fitPoints, nFitPoints)) continue;
 
                float val;
                try
                {
                    val = StsQuadraticCurvature.getCurvatureComponent(curveType);
                }
                catch(Exception e)
                {
                    StsException.systemError(this, "computeCurvature", "getCurvatureComponent failed.");
                    continue;
                }
                
                if (filterType == FILTER_ON_CHI_SQ)
                {
                	double chiSqrTest = chiSqrMultiplyer * nFitPoints;
                    double chiSqr = StsQuadraticCurvature.computeChiSquared();
                    if(chiSqr > chiSqrTest)
                    {
                        // if(StsPatchVolume.debug) StsException.systemDebug(this, "computeCurvature", "ChiSqr = " + chiSqr + " at volumeRow, volumeCol " + volumeRow + " " + volumeCol);
                        //continue;
                    	if (val > 0) val = badCurvature;
                    	if (val < 0) val = -badCurvature;
                    }
                } 
                
                curvature[patchPointRow][patchPointCol] = val;
                
                if(Math.abs(val) > curvatureTest) continue;
                // ChiSqr filtered vals not used for dataMin / dataMax & Statistics
                nValuePatchPoints++;
                sum += val;
                dataMax = Math.max(dataMax, val);
                dataMin = Math.min(dataMin, val);
            }
        }
        return nValuePatchPoints > 0;
    }

    public void drawRow(GL gl, int volumeRow, float y, float xMin, float xInc, StsColorscale colorscale, boolean is3d, boolean displayCurvature)
    {
        boolean lineStarted = false;
        float z;

        if(volumeRow < rowMin || volumeRow > rowMax) return;
        try
        {
            float x = xMin + colMin * xInc;
            int row = volumeRow - rowMin;
            gl.glDepthFunc(GL.GL_LEQUAL);
            // gl.glLineStipple(1, StsGraphicParameters.dottedLine);
            boolean displayCurvatureColor = (curvature != null && displayCurvature);
            if(!displayCurvatureColor)
                StsTraceUtilities.getPointTypeColor(patchType).setGLColor(gl);
            for (int col = 0; col < nCols; col++, x += xInc)
            {
                StsPatchPoint patchPoint = patchPoints[row][col];
                if (patchPoint != null)
                {
                    z = patchPoint.z;
                    if(displayCurvatureColor)
                    {
                        // StsTraceUtilities.getPointTypeColor(patchType).setGLColor(gl);
                        float v = curvature[row][col];
                        if(v == nullValue)
                            StsColor.BLACK.setGLColor(gl);
                        else
                            colorscale.getStsColor(colorscale.getIndexFromValue(v)).setGLColor(gl);
                    }                   
                    if (!lineStarted)
                    {
                        gl.glBegin(GL.GL_LINE_STRIP);
                        lineStarted = true;
                    }
                    if(is3d)
                        gl.glVertex3f(x, y, z);
                    else
                        gl.glVertex2f(x, z);
                    
                    if(patchPoint.rowCorrel == 0.0f && lineStarted)
                    {
                        lineStarted = false;
                        gl.glEnd();
                    }
                }
                else if (lineStarted)
                {
                    lineStarted = false;
                    gl.glEnd();
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "drawRow", e);
        }
        finally
        {
            if (lineStarted) gl.glEnd();
            // if(drawingDotted)gl.glDisable(GL.GL_LINE_STIPPLE);
        }
    }

    public void drawCol(GL gl, int volumeCol, float x, float yMin, float yInc, StsColorscale colorscale, boolean is3d, boolean displayCurvature)
    {
       boolean lineStarted = false;
       float z;

        try
        {
            if(volumeCol < colMin || volumeCol > colMax) return;
            float y = yMin + rowMin*yInc;
            int col = volumeCol - colMin;
            gl.glDepthFunc(GL.GL_LEQUAL);
            boolean displayCurvatureColor = curvature != null && displayCurvature;
            if(!displayCurvatureColor)
                StsTraceUtilities.getPointTypeColor(patchType).setGLColor(gl);
            float zMin = patchVolume.zMin;
            for (int row = 0; row < nRows; row++, y += yInc)
            {
                StsPatchPoint patchPoint = patchPoints[row][col];

                if (patchPoint != null && (z = patchPoint.z) != StsParameters.nullValue)
                {
                    if(displayCurvatureColor)
                    {
                        float v = curvature[row][col];
                        if(v == nullValue)
                            StsColor.BLACK.setGLColor(gl);
                        else
                            colorscale.getStsColor(colorscale.getIndexFromValue(v)).setGLColor(gl);
                    }
                    if (!lineStarted)
                    {
                        gl.glBegin(GL.GL_LINE_STRIP);
                        lineStarted = true;
                    }
                    if(is3d)
                        gl.glVertex3f(x, y, z);
                    else
                        gl.glVertex2f(y, z);

                    if(patchPoint.colCorrel == 0.0f && lineStarted)
                    {
                        lineStarted = false;
                        gl.glEnd();
                    }
                }
                else if (lineStarted)
                {
                    lineStarted = false;
                    gl.glEnd();
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "drawRow", e);
        }
        finally
        {
            if(lineStarted) gl.glEnd();
        }
    }

    public boolean isPatchGridNearZCursor(float z)
    {
        return z >= zMin && z <= zMax;
    }

    public void draw(GL gl, float xMin, float xInc, float yMin, float yInc, boolean displayCurvature, StsColorscale colorscale)
    {
       if(diamondStrips == null) diamondStrips = new StsDiamondStrips(this);
       if(displayCurvature && curvature != null)
       {
//           diamondStrips.setValues(curvature);
//           diamondStrips.drawSurfaceFillWithNullsAndColor(gl, colorscale);
       }
       else
       {
           StsColor patchColor = StsTraceUtilities.getPointTypeColor(patchType);
           patchColor.setGLColor(gl);
           diamondStrips.drawSurfaceFillWithNulls(gl);
       }
        // TriangleStrip.drawSurfaceFillWithNulls(gl, this, tStrips, pointsZ, tStripNormals, true);
    }

/*
    public void draw(GL gl, float xMin, float xInc, float yMin, float yInc, StsColorscale colorscale)
    {
        float xInc2 = xInc / 2.f;
        float yInc2 = yInc / 2.f;
        float y = yMin + rowMin*yInc;
        for(int row = rowMin; row <= rowMax; row++, y += yInc)
        {
            float[] rowPointsZ1 = pointsZ[row - rowMin];
            float[] rowVals1 = null;
            boolean drawColors = false;
            if (values != null)
            {
                rowVals1 = values[row - rowMin];
                drawColors = true;
            }
            else
                StsColor.GREY.setGLColor(gl);
            
            gl.glBegin(GL.GL_QUADS);
            float x = xMin + colMin*xInc;
            for (int col = colMin; col <= colMax; col++, x += xInc)
            {
                float z1 = rowPointsZ1[col - colMin];
                if (z1 == StsParameters.nullValue) continue;
                if(drawColors)
                {
                    float v = rowVals1[col - colMin];
                    StsColor color = colorscale.getStsColor(colorscale.getIndexFromValue(v));
                    color.setGLColor(gl);
                }
                gl.glVertex3f(x - xInc2, y - yInc2, z1);
                gl.glVertex3f(x - xInc2, y + yInc2, z1);
                gl.glVertex3f(x + xInc2, y + yInc2, z1);
                gl.glVertex3f(x + xInc2, y - yInc2, z1);
            }
            gl.glEnd();
        }
    }

    public void drawRowColorSurf_busted(GL gl, int row, int planeStartCol, int planeEndCol, float y, float xMin, float xInc, float yInc, float zInc, float z, float curvMin, float curvMax, StsColorscale colorscale)

    {
        int colStart = Math.max(colMin, planeStartCol);
        int colEnd = Math.min(colMax, planeEndCol);

        if (row < rowMin) return;
        if (row > rowMax - 1) return;

        if (zMin > z || zMax < z) return;

        float x = xMin + colStart * xInc;
        float[] rowPointsZ1 = pointsZ[row - rowMin];
        float[] rowPointsZ2 = pointsZ[1 + row - rowMin];
        float[] rowVals1 = null;
        float[] rowVals2 = null;
        if (values != null)
        {
            rowVals1 = values[row - rowMin];
            rowVals2 = values[1 + row - rowMin];
        }
        boolean begin = false;
        int count = 0;
        gl.glShadeModel(GL.GL_SMOOTH);
        gl.glDisable(GL.GL_CULL_FACE);
        gl.glDisable(GL.GL_POLYGON_STIPPLE);

        gl.glPointSize(1.f);
        gl.glColor4f(1.f, 1.f, 1.f, 1.f);
        if ((rowVals1 == null) || (rowVals2 == null))
        {
            return;
        }
        x = xMin + colStart * xInc;


        float v1[] = new float[3];
        float v2[] = new float[3];

        for (int col = colStart; col <= colEnd; col++, x += xInc)
        {

            float z1 = rowPointsZ1[col - colMin];
            float z2 = rowPointsZ2[col - colMin];

            {
                if ((z1 == StsParameters.nullValue) || (z2 == StsParameters.nullValue))
                {
                    //System.out.println("end ");
                    if (begin)
                    {
                        if (count == 2)
                        {
                            gl.glVertex3fv(v2, 0);
                            gl.glVertex3fv(v1, 0);
                        }

                        gl.glEnd();
                        begin = false;
                        count = 0;
                    }
                }
                else
                {
                    if (!begin)
                    {
                        //gl.glBegin(GL.GL_QUAD_STRIP);
                        gl.glBegin(GL.GL_LINES);
                        begin = true;
                        count = 0;
                    }
                    if (z2 == StsParameters.nullValue)
                    {
                        v1[0] = v2[0] = x;
                        v1[1] = v2[1] = y;
                        v1[2] = v2[2] = z1;
                        gl.glVertex3fv(v2, 0);
                        gl.glVertex3fv(v1, 0);
                        count += 2;


                    }
                    else if (z1 == StsParameters.nullValue)
                    {
                        v1[0] = v2[0] = x;
                        v1[1] = v2[1] = y + yInc;
                        v1[2] = v2[2] = z2;
                        gl.glVertex3fv(v2, 0);
                        gl.glVertex3fv(v1, 0);
                        count += 2;


                    }
                    else
                    {

                        v2[0] = x;
                        v2[1] = y + yInc;
                        v2[2] = z2;
                        gl.glVertex3fv(v2, 0);

                        v1[0] = x;
                        v1[1] = y;
                        v1[2] = z1;
                        gl.glVertex3fv(v1, 0);

                        count += 2;
                        //System.out.println(" both "+count);
                        //System.out.println("    "+v1[0]+" "+v1[1]+" "+v1[2]);
                        //System.out.println("    "+v2[0]+" "+v2[1]+" "+v2[2]);
                    }
                }

            }
            if (begin)
            {
                if (count == 2)
                {

                    gl.glVertex3fv(v2, 0);
                    gl.glVertex3fv(v1, 0);
                }

                gl.glEnd();
            }
        }
    }

    transient boolean ppChecked = false;
    transient boolean hasPP = false;

    boolean hasPointParams(GL gl)
    {
        if (ppChecked) return hasPP;

        if (gl.isExtensionAvailable("GL_ARB_point_parameters"))
            hasPP = true;

        if (gl.isExtensionAvailable("GL_EXT_point_parameters"))
            hasPP = true;

        ppChecked = true;

        return hasPP;
    }

    double sqr(double a) { return (a * a); }

    private void setPointParams(GL gl)
    {
        if (!hasPointParams(gl)) return;
        int error = gl.glGetError();
        float maxSize = 50.f;

        gl.glPointParameterf(GL.GL_POINT_SIZE_MIN, 2.0f);
        double[] projectionMatrix = new double[16];
        int[] viewport = new int[4];
        gl.glGetDoublev(GL.GL_PROJECTION_MATRIX, projectionMatrix, 0);
        gl.glGetIntegerv(GL.GL_VIEWPORT, viewport, 0);
        double H = viewport[2];
        double h = 2.0 / projectionMatrix[0];
        double D0 = Math.sqrt(2.0 * H / h);
        double k = 1.0 / (1.0 + 2 * sqr(1. / projectionMatrix[0]));
        float[] atten = new float[3];
        atten[0] = 1.f;
        atten[1] = 0.0f;
        k /= 500.f;
        atten[2] = (float) sqr(1 / D0) * (float) k;
        //System.out.println(" atten "+atten[2]);
        //atten[2] = 0.000001f;

        gl.glPointParameterfv(GL.GL_DISTANCE_ATTENUATION_EXT, atten, 0);
        error = gl.glGetError();
        if (error != 0)
        {
            GLU glu = new GLU();
            System.out.println("pointParams err code " + error + " " + glu.gluErrorString(error));
        }
        System.out.println("atten " + atten[2]);
        float[] ft = new float[1];
        gl.glGetFloatv(GL.GL_POINT_SIZE_MAX, ft, 0);
        System.out.println("point max " + ft[0]);
        gl.glPointSize(ft[0] > maxSize ? maxSize : ft[0]);
        gl.glPointParameterf(GL.GL_POINT_SIZE_MAX, ft[0] > maxSize ? maxSize : ft[0]);
    }
*/
    public void drawRowVox(GL gl, float yMin, float yInc, float xMin, float xInc, StsColorscale colorscale)
    {
        float y = yMin + rowMin*yInc;
        for(int row = rowMin; row <= rowMax; row++, y += yInc)
        {
            float[] rowPointsZ = pointsZ[row - rowMin];

            if (curvature == null) continue;
            float[] rowVals = curvature[row - rowMin];
            if (rowVals == null) continue;
            gl.glPointSize(3.f);
            gl.glBegin(GL.GL_POINTS);
            float x = xMin + colMin*xInc;
            for (int col = colMin; col <= colMax; col++, x += xInc)
            {
                float z = rowPointsZ[col - colMin];
                Color color = colorscale.getColor(colorscale.getIndexFromValue(rowVals[col - colMin]));
                float colorsf[] = new float[4];
                color.getRGBComponents(colorsf);
                gl.glColor4fv(colorsf, 0);
                if (z != StsParameters.nullValue)
                    gl.glVertex3f(x, y, z);
            }
            gl.glEnd();
        }
    }

    public int size()
    {
        return (rowMax - rowMin + 1) * (colMax - colMin + 1);
    }

    public String toString()
    {
        String rowColString = super.toString();
        if (connectionsList != null)
        {
            if(connectionsList.size() > 0)
            {
                StsPatchConnection connection = connectionsList.get(0);
                StsPatchPoint patchPoint = connection.newPoint;
                return "id: " + id + " nConnections " + connectionsList.size() + " " + rowColString + " zFirst: " + patchPoint.value;
            }
            else
                return "id: " + id + " " + rowColString + " no points";
        }
        else if (patchPoints != null)
            return "id: " + id + " idFinal: " + idFinal + " nPatchPoints " + nPatchPoints + " " + rowColString + " zMin: " + zMin + " zMax: " + zMax;
        else
            return "id: " + id + " idFinal: " + idFinal + " nPatchPoints " + nPatchPoints + " " +  rowColString;
    }

    public int fillHistogram(float[] data, int nValues)
    {
        for(int row = 0; row < nRows; row++)
        {
            for(int col = 0; col < nCols; col++)
            {
                float value = curvature[row][col];
                if(value != nullValue)
                    data[nValues++] = value;
                if(nValues == data.length)
                    return nValues;
            }
        }
        return nValues;
    }

    public float getXMin() { return patchVolume.xMin; }
    public float getXMax() { return patchVolume.xMax; }
    public float getYMin() { return patchVolume.yMin; }
    public float getYMax() { return patchVolume.yMax; }
    public float getXInc() { return patchVolume.xInc; }
    public float getYInc() { return patchVolume.yInc; }
    public float getRowCoor(float[] xyz) { return patchVolume.getRowCoor(xyz); }
    public float getColCoor(float[] xyz) { return patchVolume.getRowCoor(xyz); }
    public double getXOrigin() { return patchVolume.xOrigin; }
    public double getYOrigin() { return patchVolume.yOrigin; }
    public float getXSize() { return patchVolume.getXSize(); }
    public float getYSize() { return patchVolume.getYSize(); }
    public float getAngle() { return patchVolume.getAngle(); }
    public float getXCoor(float rowF, float colF) { return patchVolume.getXCoor(rowF, colF); }
    public float getYCoor(float rowF, float colF) { return patchVolume.getYCoor(rowF, colF); }
    public StsPoint getPoint(int row, int col)
    {
        float[] xyz = getXYZorT(row, col);
        return new StsPoint(xyz);
    }

    public float[] getXYZorT(int row, int col)
    {
        float[] xy = patchVolume.getXYCoors(row, col);
        StsPatchPoint patchPoint = getPatchPoint(row, col);
        float z = patchPoint.z;
        return new float[] { xy[0], xy[1], z };
    }

    public StsPoint getPoint(float rowF, float colF) { return null; } // not used
    public float[] getXYZorT(float rowF, float colF) { return null; } // not used
    public float[] getNormal(int row, int col) { return null; } // not used
    public float[] getNormal(float rowF, float colF) { return null; } // not used
    public void checkConstructGridNormals() { return; } // not used
    public float getZInc() { return 0.0f; } // not used

    public float getZMin() { return zMin; }
    public float getZMax() { return zMax; }
    public String getLabel() { return toString(); }
    public float interpolateBilinearZ(StsPoint point, boolean computeIfNull, boolean setPoint) { return 0.0f; } // not used:  yet
    public float interpolateBilinearZ(StsGridPoint gridPoint, boolean computeIfNull, boolean setPoint) { return 0.0f; } // not used: yet

    public float getComputePointZ(int row, int col) { return pointsZ[row][col]; }

    public boolean toggleSurfacePickingOn() { return false; }
    public void toggleSurfacePickingOff() { }
    public String getName() { return "patchGrid-" + idFinal; }
    public StsGridPoint getSurfacePosition(StsMouse mouse, boolean display, StsGLPanel3d glPanel3d) { return null; }
    public void setIsVisible(boolean isVisible) { this.isVisible = isVisible; }
    public boolean getIsVisible() { return isVisible; }

    public boolean hasRowLink(int row, int col)
    {
        return getRowCorrel(row, col) > patchVolume.minLinkCorrel;
    }
    public boolean hasColLink(int row, int col)
    {
        return getColCorrel(row, col) > patchVolume.minLinkCorrel;
    }

    public StsPatchPoint getPatchPoint(int row, int col)
    {
        if(patchPoints == null) return null;
        return patchPoints[row-rowMin][col-colMin];
    }

    public StsPatchPoint[] getRowPatchPoints(int volumeRow)
    {
        if(patchPoints == null) return null;
        return patchPoints[volumeRow-rowMin];
    }

    public float getPointZ(int volumeRow, int volumeCol)
    {
        if(!isInsideRowCol(volumeRow, volumeCol)) return nullValue;
        StsPatchPoint patchPoint = patchPoints[volumeRow-rowMin][volumeCol - colMin];
        if(patchPoint == null) return nullValue;
        return patchPoint.z;
    }

    public float getRowCorrel(int patchRow, int patchCol)
    {
        StsPatchPoint patchPoint = patchPoints[patchRow][patchCol];
        if(patchPoint == null) return 0.0f;
        return patchPoint.rowCorrel;
    }

    public float getColCorrel(int patchRow, int patchCol)
    {
        StsPatchPoint patchPoint = patchPoints[patchRow][patchCol];
        if(patchPoint == null) return 0.0f;
        return patchPoint.colCorrel;
    }

    public float getCurvature(int volumeRow, int volumeCol, float dataMin, float dataMax)
    {
        if(!isInsideRowCol(volumeRow, volumeCol))
        {
            StsException.systemError(this, "getValue", "volumeRow or volumeCol not inside patch");
            return 0.0f;
        }
        if(curvature == null) return 0.0f;
        float value = curvature[volumeRow - rowMin][volumeCol - colMin];
        if(value == nullValue) return 0.0f;
        if (value < dataMin) return dataMin;
        else if(value > dataMax) return dataMax;
        else return value;
    }
}
