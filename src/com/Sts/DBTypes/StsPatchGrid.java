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
    public int id;
     /** the final id of this patch in gridList; this patch may be a layer from an original patch whose id is id */
    int idFinal;
    /** type of this grid: Min, Max +Zero-crossing, or -Zero-crossing */
    byte patchType;
    boolean isVisible = true;
    StsPatchVolume patchVolume;
    int nPatchPoints;
    float dataMin;
    float dataMax;
    float zMin = StsParameters.largeFloat;
    float zMax = -StsParameters.largeFloat;

    int nRows = 0;
    int nCols = 0;

    float[][] pointsZ;
    float[][] rowCorrels;
    float[][] colCorrels;

    transient StsPatchPoint[][] patchPoints = null;
    /** flag indicating this patchGrid has been added to current patchVolume.rowGrids hashmap list; avoids overhead of trying to re-add */
    transient public boolean rowGridAdded = false;

    transient RowColGrid rowColGrid = new RowColGrid(rowMin, rowMax, colMin, colMax);

    transient float[][] value;

    transient float[][] curvature;

    transient StsDiamondStrips diamondStrips;

    transient public int nValuePatchPoints = 0;
    transient public double sum = 0;
    /** points on this grid are in a hashMap with a hash key whose value is the grid index: col + row*patchVolume.nCols */
    // transient HashMap<Integer, StsPatchPoint> patchPointsHashMap;
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

    static final int largeInt = Integer.MAX_VALUE;

    static final int NO_PATCH_DEBUG = -1;
    /** various debug print of patches in rowGrid and prevRowGrid arrays; check if this is not -1 and prints if id matches this value.  Make this -1 if you want no debug prints. */
    static final int debugPatchInitialID = NO_PATCH_DEBUG;
    /** debugPatchID may change if original patch is merged into a new one; when this occurs, set debugCurrentPatchID to new one and track it */
    static int debugPatchID = debugPatchInitialID;

    static final int debugPointRow = 8;
    static final int debugPointCol = 8;

    static final boolean debug = debugPatchInitialID != -1;
    static final boolean debugPoint = debugPointRow != -1 && debugPointCol != -1;

    public StsPatchGrid()
    {
    }
    
    private StsPatchGrid(StsPatchVolume patchVolume, byte patchType)
    {
        this.patchVolume = patchVolume;
        this.patchType = patchType;
    }

    static public StsPatchGrid construct(StsPatchVolume patchVolume, byte patchType)
    {
        StsPatchGrid patchGrid = new StsPatchGrid(patchVolume, patchType);
        patchGrid.setup();
        return patchGrid;
    }

    static public void staticInitialize()
    {
        nextPatchID = 0;
        nextFinalPatchID = 0;
        debugPatchID = debugPatchInitialID;
    }

    private void setup()
    {
        id = nextPatchID++;
        if(debugPatchID != -1 && id == debugPatchID)
            //if(debug())
            StsException.systemDebug(this, "setup", "patch " + id + " initialized");
    }

    public boolean debug()
    {
        return debugPatchID != -1 && id == debugPatchID;
    }

    /** give two different size grids with this intersection, check if we can merge the two;
     *  possible only if they don't have common occupancy of any row-col location
     * @param patchGrid1 first of patchGrids to merge
     * @param patchGrid2 second of patchGrids to merge
     */
    static public boolean mergePatchPointsOK(StsPatchGrid patchGrid1, StsPatchGrid patchGrid2)
    {
        // check for any overlap between this grid and patchPointGrid
        RowColGrid intersectGrid = patchGrid1.rowColGridIntersect(patchGrid2);
        int r1 = intersectGrid.rowMin - patchGrid1.rowMin;
        int r2 = intersectGrid.rowMin - patchGrid2.rowMin;
        int c1 = intersectGrid.colMin - patchGrid1.colMin;
        int c2 = intersectGrid.colMin - patchGrid2.colMin;
        for(int row = 0; row < intersectGrid.nRows; row++, r1++, r2++)
        {
            int cc1 = c1;
            int cc2 = c2;
            for (int col = 0; col < intersectGrid.nCols; col++, cc1++, cc2++)
            {
                if(patchGrid1.patchPoints[r1][cc1] != null && patchGrid2.patchPoints[r2][cc2] != null)
                    return false;
            }
        }
        return true;
    }

    /** We wish to merge the points from removedGrid into this one.  Copy both sets of points to a new grid which is union of two.
     *  reset the removedGrid patchPoints.id to this id
     * @param removedGrid  newPatchGrid to be merged to this (otherPatchGrid).
     * @return true if merged successfully
     */
    boolean mergePatchPoints(StsPatchGrid removedGrid)
    {
        if(debug && (id == debugPatchID || removedGrid.id == debugPatchID))
            StsException.systemDebug(this, "mergePatchPoints", "merging patch id: " + removedGrid.id + " to patch id: " + id);

        RowColGrid union = rowColGridUnion(removedGrid); //union of this and newPatchGrid
        StsPatchPoint[][] newMergedPatchPoints = new StsPatchPoint[union.nRows][union.nCols]; // create empty mergedPoints grid
        if(!copyPatchPointsTo(newMergedPatchPoints, union)) return false; // copy this patchPoints to mergedPoints
        if(!removedGrid.copyPatchPointsTo(newMergedPatchPoints, union)) return false;
        removedGrid.resetPatchPointsGrid(this);
        resetPatchPoints(union, newMergedPatchPoints);
        nPatchPoints += removedGrid.nPatchPoints;
        if(debugPatchID != -1 && removedGrid.id == debugPatchID)
            debugPatchID = id;
        return true;
    }

    /** this patchGrid is being merged to another patchGrid with this ID. Reset all the merged patchPoints to this id.
     *
     * @param newPatchGrid patch which this patch is being merged into
     */
    private void resetPatchPointsGrid(StsPatchGrid newPatchGrid)
    {
        for(int row = 0; row < nRows; row++)
            for(int col = 0; col < nCols; col++)
                if(patchPoints[row][col] != null) patchPoints[row][col].setPatchGrid(newPatchGrid);
    }

    RowColGrid rowColGridIntersect(StsPatchGrid otherGrid)
    {
        int rowMin, rowMax, colMin, colMax;
        rowMin = Math.max(this.rowMin, otherGrid.rowMin);
        rowMax = Math.min(this.rowMax, otherGrid.rowMax);
        colMin = Math.max(this.colMin, otherGrid.colMin);
        colMax = Math.min(this.colMax, otherGrid.colMax);
        return new RowColGrid(rowMin, rowMax, colMin, colMax);
    }

    RowColGrid rowColGridUnion(StsPatchGrid otherGrid)
    {
        int rowMin, rowMax, colMin, colMax;
        rowMin = Math.min(this.rowMin, otherGrid.rowMin);
        rowMax = Math.max(this.rowMax, otherGrid.rowMax);
        colMin = Math.min(this.colMin, otherGrid.colMin);
        colMax = Math.max(this.colMax, otherGrid.colMax);
        return new RowColGrid(rowMin, rowMax, colMin, colMax);
    }

    private RowColGrid getRowColGrid()
    {
        return new RowColGrid(rowMin, rowMax, colMin, colMax);
    }

    class RowColGrid
    {
        int rowMin = largeInt;
        int rowMax = -largeInt;
        int colMin = largeInt;
        int colMax = -largeInt;
        int nRows = 0, nCols = 0;

        RowColGrid(int rowMin, int rowMax, int colMin, int colMax)
        {
            this.rowMin = rowMin;
            this.rowMax = rowMax;
            this.colMin = colMin;
            this.colMax = colMax;
            nRows = rowMax - rowMin + 1;
            nCols = colMax - colMin + 1;
        }

        public String toString()
        {
            return new String("rowMin: " + rowMin + " rowMax: " + rowMax + " colMin: " + colMin + " colMax: " + colMax);
        }
    }

    boolean patchPointOverlaps(StsPatchPoint patchPoint)
    {
        try
        {
            if (!contains(patchPoint)) return false;
            return getPatchPoint(patchPoint.row, patchPoint.col) != null;

        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "patchPointOverlaps", "Failed for point: " + patchPoint.toString() +
                    " on grid: " + toString(), e);
            return false;
        }
    }

    void addPatchPoint(StsPatchPoint patchPoint)
    {
        if(debug && id == debugPatchID)
        //       if(debug && id == debugPatchID && debugPointRow == patchPoint.row && debugPointCol == patchPoint.col)
            StsException.systemDebug(this, "addPatchPoint", " patchPoint " + patchPoint.toString());

        if(patchPoints == null)
            initializePatchPoints(patchPoint);
        else
            checkAdjustGrid(patchPoint);
        if(!contains(patchPoint))
        {
            StsException.systemError(this, "addPatchPoint", "pointGrid " + rowColGrid.toString() + " doesn't contain point " + patchPoint.toString());
            return;
        }
        patchPoints[patchPoint.row - rowMin][patchPoint.col - colMin] = patchPoint;
        patchPoint.setPatchGrid(this);
        nPatchPoints++;
    }

    void checkAdjustGrid(StsPatchPoint patchPoint)
    {
        int row = patchPoint.row;
        int col = patchPoint.col;

        int rowMinNew = rowMin, rowMaxNew = rowMax, colMinNew = colMin, colMaxNew = colMax;

        boolean gridChanged = false;
        if (row < this.rowMin)
        {
            rowMinNew = row;
            gridChanged = true;
        }
        if (row > this.rowMax)
        {
            rowMaxNew = row;
            gridChanged = true;
        }
        if (col < this.colMin)
        {
            colMinNew = col;
            gridChanged = true;
        }
        if (col > this.colMax)
        {
            colMaxNew = col;
            gridChanged = true;
        }

        if (!gridChanged) return;

        RowColGrid newRowColGrid = new RowColGrid(rowMinNew, rowMaxNew, colMinNew, colMaxNew);
        copyResetRowColGrid(newRowColGrid);
    }

    void copyResetRowColGrid(RowColGrid newRowColGrid)
    {
        if(debug && debugPatchID != -1 && id == debugPatchID)
            StsException.systemDebug(this, "copyResetRowColGrid", "grid reset from " + rowColGrid + " to " + newRowColGrid);
        StsPatchPoint[][] newPatchPoints = copyPatchPoints(newRowColGrid);
        resetPatchPoints(newRowColGrid, newPatchPoints);
    }

    StsPatchPoint[][] copyPatchPoints(RowColGrid newRowColGrid)
    {
        if (patchPoints == null) return null;
        StsPatchPoint[][] newPatchPoints = new StsPatchPoint[newRowColGrid.nRows][newRowColGrid.nCols];
        if(!copyPatchPointsTo(newPatchPoints, newRowColGrid))
            return null;
        else
            return newPatchPoints;
    }

    boolean copyPatchPointsTo(StsPatchPoint[][] newPatchPoints, RowColGrid newRowColGrid)
    {
        int row = -1, newRow = -1;
        int col = -1, newCol = -1;
        try
        {
            int rowStart = rowMin - newRowColGrid.rowMin;
            int colStart = colMin - newRowColGrid.colMin;
            for (row = 0, newRow = rowStart; row < nRows; row++, newRow++)
                for(col = 0, newCol = colStart; col < nCols; col++, newCol++)
                {
                    if(patchPoints[row][col] != null)
                    {
                        if(newPatchPoints[newRow][newCol] != null)
                            StsException.systemError(this, "copyPatchPointsTo", "Failed copying patch " + rowColGrid.toString() +
                                    " row: " + row + " col: " + col + " to new grid " + newRowColGrid.toString() +
                                     "row: " + newRow + "col: " + newCol);
                        else
                            newPatchPoints[newRow][newCol] = patchPoints[row][col];
                    }
                }
            return true;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "copyPatchPointsTo", "Failed copying patch " + rowColGrid.toString() + " row: " + row +
                " to new grid " + newRowColGrid.toString() + "row: " + newRow, e);
            return false;
        }
    }

    void resetPatchPoints(RowColGrid newRowColGrid, StsPatchPoint[][] newPatchPoints)
    {
        initializeRowColGrid(newRowColGrid);
        patchPoints = newPatchPoints;
    }

    boolean contains(RowColGrid newRowColGrid)
    {
        return newRowColGrid.rowMin >= rowMin && newRowColGrid.rowMax <= rowMax &&
                newRowColGrid.colMin >= colMin && newRowColGrid.colMax <= colMax;
    }

    boolean contains(StsPatchPoint patchPoint)
    {
        int row = patchPoint.row;
        int col = patchPoint.col;
        return row >= rowMin && row <= rowMax && col >= colMin && col <= colMax;
    }

    private void initializePatchPoints(StsPatchPoint patchPoint)
    {
        initializeRowColGrid(patchPoint);
        patchPoints = new StsPatchPoint[nRows][nCols];
    }

    /** Called to add a point which doesn't belong to a grid to this patchGrid as long as it doesn't overlap.
     *  If it does overlap, create a new patchGrid and add this point and clone of connectedPatchPoint which
     *  belongs to this grid.
     */
    public StsPatchGrid checkAddPatchPoint(StsPatchPoint patchPoint, StsPatchPoint connectedPatchPoint)
    {
        if(patchPointOverlaps(patchPoint))
        {
            StsPatchGrid newPatchGrid = StsPatchGrid.construct(patchVolume, patchType);
            newPatchGrid.addPatchPoint(patchPoint);
            newPatchGrid.addPatchPoint(connectedPatchPoint.clone());
            return newPatchGrid;
        }
        else
        {
            addPatchPoint(patchPoint);
            return this;
        }
    }

    public final void addCorrelation(StsPatchPoint otherPatchPoint, StsPatchPoint newPatchPoint, float correl)
    {
        if (otherPatchPoint.row == newPatchPoint.row)
            addRowCorrelation(otherPatchPoint,  newPatchPoint, correl);
        else
            addColCorrelation(otherPatchPoint,  newPatchPoint, correl);
    }

    public final void addRowCorrelation(StsPatchPoint otherPatchPoint, StsPatchPoint newPatchPoint, float correl)
    {
        if(debug)
        {
            if(otherPatchPoint.z == nullValue)
                StsException.systemError(this, "addRowCorrelation", "Has row correl at row: " + otherPatchPoint.row +
                " col: " + otherPatchPoint.col + " but first point is null.");
            if(newPatchPoint.z == nullValue)
                StsException.systemError(this, "addRowCorrelation", "Has row correl at row: " + otherPatchPoint.row +
                        " col: " + otherPatchPoint.col + " but second point is null.");
        }
        otherPatchPoint.rowCorrel = correl;
    }

    public final void addColCorrelation(StsPatchPoint otherPatchPoint, StsPatchPoint newPatchPoint, float correl)
    {
        if(debug)
        {
            if(otherPatchPoint.z == nullValue)
                StsException.systemError(this, "addColCorrelation", "Has col correl at row: " + otherPatchPoint.row +
                        " col: " + otherPatchPoint.col + " but first point is null.");
            if(newPatchPoint.z == nullValue)
                StsException.systemError(this, "addColCorrelation", "Has col correl at row: " + otherPatchPoint.row +
                        " col: " + otherPatchPoint.col + " but second point is null.");
        }
        otherPatchPoint.colCorrel = correl;
    }

    public void clear()
    {
        initializeRowColGrid();
        patchPoints = null;
        if(debugPatchID != -1 && id == debugPatchID)
        //if(debug())
            StsException.systemDebug(this, "clear", "clearing patch: " + toString());
    }

    private void initializeRowColGrid()
    {
        initializeRowColGrid(largeInt, -largeInt, largeInt, -largeInt);
    }

    private void initializeRowColGrid(int rowMin, int rowMax, int colMin, int colMax)
    {
        initializeRowColGrid(new RowColGrid(rowMin, rowMax, colMin, colMax));
    }

    private void initializeRowColGrid(StsPatchPoint patchPoint)
    {
        int row = patchPoint.row;
        int col = patchPoint.col;
        initializeRowColGrid(row, row, col, col);
    }

    private void initializeRowColGrid(RowColGrid newRowColGrid)
    {
        if(debugPatchID != -1 && id == debugPatchID)
            //if(debug())
            StsException.systemDebug(this, "initializeRowColGrid", "Reset rowColGrid from " + rowColGrid + " to " + newRowColGrid);
        rowMin = newRowColGrid.rowMin;
        rowMax = newRowColGrid.rowMax;
        colMin = newRowColGrid.colMin;
        colMax = newRowColGrid.colMax;
        nRows = rowMax - rowMin + 1;
        nCols = colMax - colMin + 1;
        this.rowColGrid = newRowColGrid;
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

    public void resetIndex(int index)
    {
        if(debugPatchID != -1 && id == debugPatchID)
            //if(debug())
            StsException.systemDebug(this, "resetIndex", "debugPatch final ID being reset from " + idFinal + " to " + index);
        idFinal = index;
    /*
        for(int row = 0; row < nRows; row++)
            for(int col = 0; col < nCols; col++)
                if(patchPoints[row][col] != null)
                    patchPoints[row][col].patchID = idFinal;
    */
    }

    public void finish()
    {
        if(debug && id == debugPatchID)
            StsException.systemDebug(this, "finish", "for patch " + toString());
        pointsZ = new float[nRows][nCols];
        rowCorrels = new float[nRows][nCols];
        colCorrels = new float[nRows][nCols];
        for (int row = 0; row < nRows; row++)
        {
            for (int col = 0; col < nCols; col++)
            {
                if (patchPoints[row][col] != null)
                {
                    float z = patchPoints[row][col].z;
                    pointsZ[row][col] = z;
                    zMin = Math.min(zMin, z);
                    zMax = Math.max(zMax, z);
                    rowCorrels[row][col] = patchPoints[row][col].rowCorrel;
                    colCorrels[row][col] = patchPoints[row][col].colCorrel;
                }
                else
                {
                    pointsZ[row][col] = nullValue;
                    rowCorrels[row][col] = 0.0f;
                    colCorrels[row][col] = 0.0f;
                }
            }
        }
        if(!StsPatchVolume.debug) patchPoints = null;
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

    public float[][] getPointsZ() { return pointsZ; }

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
                if (pointsZ[row][col] != nullValue)
                {
                    colStart = col;
                    break;
                }
            }
            for (int col = nCols - 1; col > 0; col--)
            {
                if (pointsZ[row][col] != nullValue)
                {
                    colEnd = col;
                    break;
                }
            }

            for(int col = colStart; col <= colEnd; col++)
                if(pointsZ[row][col] != nullValue) nActualUsed++;

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
        int col = -1;

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
            for (col = 0; col < nCols; col++, x += xInc)
            {
                z = pointsZ[row][col];
                if (z != nullValue)
                {
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
                    
                    if(rowCorrels[row][col] == 0.0f && lineStarted)
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
            StsException.outputWarningException(this, "drawRow", "Failed for patchGrid " + id + " at row: " + volumeRow + " col: " + col, e);
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
       int row = -1;

        try
        {
            if(volumeCol < colMin || volumeCol > colMax) return;
            float y = yMin + rowMin*yInc;
            int col = volumeCol - colMin;
            gl.glDepthFunc(GL.GL_LEQUAL);
            boolean displayCurvatureColor = curvature != null && displayCurvature;
            if(!displayCurvatureColor)
                StsTraceUtilities.getPointTypeColor(patchType).setGLColor(gl);
            for (row = 0; row < nRows; row++, y += yInc)
            {
                z = getPointZ(row, col);

                if (z != StsParameters.nullValue)
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

                    if(colCorrels[row][col] == 0.0f && lineStarted)
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
            StsException.outputWarningException(this, "drawCol", "Failed for patchGrid " + id + " at row: " + row + " col: " + volumeCol, e);
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
        if (patchPoints != null)
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
    public StsPoint getPoint(int volumeRow, int volumeCol)
    {
        float[] xyz = getXYZorT(volumeRow, volumeCol);
        return new StsPoint(xyz);
    }

    public float[] getXYZorT(int volumeRow, int volumeCol)
    {
        float[] xy = patchVolume.getXYCoors(volumeRow, volumeCol);

        float z = pointsZ[volumeRow - rowMin][volumeCol - colMin];
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

    public float getPointZ(int patchRow, int patchCol)
    {
        if(pointsZ == null)
            return nullValue;
        if(!isInsidePatchRowCol(patchRow, patchCol))
            return nullValue;
        return pointsZ[patchRow][patchCol];
    }

    public boolean isInsidePatchRowCol(int row, int col)
    {
        return row >= 0 && row < nRows && col >= 0 && col < nCols;
    }

    public float getRowCorrel(int patchRow, int patchCol)
    {
        if(rowCorrels == null)
            return 0.0f;
        if(!isInsidePatchRowCol(patchRow, patchCol))
            return 0.0f;
        return rowCorrels[patchRow][patchCol];
    }

    public float getColCorrel(int patchRow, int patchCol)
    {
        if(colCorrels == null)
            return 0.0f;
        if(!isInsidePatchRowCol(patchRow, patchCol))
            return 0.0f;
        return colCorrels[patchRow][patchCol];
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
