package com.Sts.Utilities.Interpolation;

import com.Sts.DBTypes.*;
import com.Sts.Reflect.*;
import com.Sts.Utilities.*;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jan 30, 2007
 * Time: 9:02:29 AM
 * To change this template use File | Settings | File Templates.
 */
abstract public class StsRadialInterpolation
{
    public StsPreStackLineSet lineSet;

    /** A regular grid volume will be nRows x nCols. So we set nRows and nCols;
     *  A volume which is a collection of lines will have a row per line, and nCols per row. So we set nRows and nColsInRow array.
     */
    public int nRows = 0;
    public int nCols = 0;
    public int[] nColsInRow;
    public int nTotalPoints = 0;
    DistanceSq[][][] gridWeights;
    public int[][] nGridWeights = null;
    public boolean[][] isNeighbor; // true: is in the neighborhood of a point

    public boolean[][] gridChanged; // true: weights at this point have changed
    public boolean[] rowChanged; // point in row has changed
    public boolean[] colChanged; // point in col has changed
    public boolean changed;

    int nMinPoints;
    int nMaxPoints;//	int maxWeights;
    float maxDistSq;
    float neighborRadius;
    public boolean upToDate = false;
    int nUpdatedPoints;
    StsMethod updateMethod;
    boolean modify = true;
    int currentMinNWeights = StsParameters.largeInt;

    ArrayList[] dataPoints;
    ArrayList[] newDataPoints;

    static final boolean debug = false;
 
    abstract public void initialize();
    abstract public void run();
    public abstract void runUpdate(Object object, int row, int col);

    static public StsRadialInterpolation constructor(StsPreStackLineSet lineSet, int nMinPoints, int nMaxPoints, float neighborRadius)
    {
        if(lineSet instanceof StsPreStackLineSet3d)
            return new StsSpiralRadialInterpolation((StsPreStackLineSet3d)lineSet, nMinPoints, nMaxPoints, neighborRadius);
        else
            return new StsLinesRadialInterpolation((StsPreStackLineSet2d)lineSet, nMinPoints, nMaxPoints, neighborRadius);
    }

    public void setUpdateMethod(StsMethod updateMethod)
    {
        this.updateMethod = updateMethod;
    }

    protected boolean addWeight(DataPoint dataPoint, int row, int col, double distSq)
    {
//		if(mainDebug && row == debugRow && col == debugCol)
//		{
//			System.out.println("break");
//		}
//		if(distSq > maxDistSq) return false;

        if (nGridWeights[row][col] == -1) return false; // this is a data point
        DistanceSq[] pointWeights = gridWeights[row][col];
        int nWeights = nGridWeights[row][col];
        DistanceSq pointWeight = new DistanceSq(row, col, dataPoint, distSq);
        if (nWeights == 0)
        {
            pointWeights = new DistanceSq[nMaxPoints];
            pointWeights[0] = pointWeight;
            nGridWeights[row][col] = 1;
            gridWeights[row][col] = pointWeights;
            setGridStatus(row, col, distSq);
            currentMinNWeights = Math.min(currentMinNWeights, 1);
            return true;
        }
        else if (nWeights > nMaxPoints)
            return false;
        else
        {
            for (int n = nWeights - 1; n >= 0; n--)
            {
                if (distSq >= pointWeights[n].distanceSq)
                    return insert(pointWeight, pointWeights, nWeights, row, col, n + 1);
                else if (n == 0)
                    return insert(pointWeight, pointWeights, nWeights, row, col, 0);
            }
            return false;
        }
    }

    /** When a new point is added, it is added to the newDataPoint list and upToDate is set to false;
     *  If not upToDate, we run all new points in newDataPoint list; otherwise we runUpdate which flags any
     *  row-col locations changed by this added object.
     * @param object recent object added
     * @param row row location of this object
     * @param col col location of this object
     */
    public void updateDataPoint(Object object, int row, int col)
    {
//        setGridStatus(row, col, 0.0f);
        if(!upToDate)
            run();
        else
            runUpdate(object, row, col);
//		run();
    }

    protected void setGridStatus(int row, int col, double distSq)
    {
        if (distSq <= maxDistSq) isNeighbor[row][col] = true;
        gridChanged[row][col] = true;
        rowChanged[row] = true;
        if(colChanged != null) colChanged[col] = true;
        changed = true;

       if(updateMethod == null || nGridWeights[row][col] < 0) return; // if this is a dataPoint, don't update
        try
        {
            updateMethod.invokeInstanceMethod(new Object[]{new Integer(row), new Integer(col)});
        }
        catch (Exception e)
        {
            StsException.outputException("StsSpiralRadialInterpolation.setGridStatus() failed invoking updateMethod for row " + row + " col " + col,
                    e, StsException.WARNING);
        }
    }

    /**
     * We have updated this row-col on this row. If the col-line it crosses has not
     * been changed (or was updated and flagged as not changed), then we can clear
     * this row-col location by setting gridChanged to false.
     */
    public void setRowStatus(int row, int col)
    {
        if(colChanged != null)
            if (!colChanged[col]) gridChanged[row][col] = false;
        else
            gridChanged[row][col] = false;
    }

    /**
     * We have updated this row-col on this col. If the row-line it crosses has not
     * been changed (or was updated and flagged as not changed), then we can clear
     * this row-col location by setting gridChanged to false.
     */
    public void setColStatus(int row, int col)
    {
        if (!rowChanged[row]) gridChanged[row][col] = false;
    }

    public void setRowStatus(int row)
    {
        rowChanged[row] = false;
    }

    public void setColStatus(int col)
    {
        if(colChanged == null) return;
        colChanged[col] = false;
    }

    private boolean insert(DistanceSq pointWeight, DistanceSq[] pointWeights, int nWeights, int row, int col, int index)
    {
        if (index >= nMaxPoints) return false;
        // move points up if there is room in list
        if (nWeights < nMaxPoints)
        {
            nWeights++;
            currentMinNWeights = Math.min(currentMinNWeights, nWeights);
            nGridWeights[row][col] = nWeights;
        }
        for (int i = nWeights - 1; i > index; i--)
            pointWeights[i] = pointWeights[i - 1];
        pointWeights[index] = pointWeight;
        setGridStatus(row, col, pointWeight.distanceSq);
        return true;
    }

    /*
        public boolean isPlaneOK(int dirNo, int nPlane)
            {
                if(dirNo == StsParameters.XDIR)
                    return colOK[nPlane];
                else if(dirNo == StsParameters.YDIR)
                    return rowOK[nPlane];
                else
                    return true;
            }
        */
    public Weights getWeights(int row, int col)
    {
        if (row == -1 || col == -1) return null;
        if (nGridWeights == null) return null;
        int nWeights = nGridWeights[row][col];
        if (nWeights >= 0 && nWeights < nMinPoints) return null;
        DistanceSq[] distanceSqs = this.gridWeights[row][col];
        return new Weights(nWeights, distanceSqs);
    }

    public boolean isDataPoint(int row, int col)
    {
        return nGridWeights[row][col] == -1;
    }

    public int  getNGridWeights(int row, int col)
    {
        return nGridWeights[row][col];
    }

    public boolean rowHasDataPoint(int row)
    {
        int nCols = nGridWeights[row].length;
        for(int col = 0; col < nCols; col++)
            if(nGridWeights[row][col] == -1) return true;
        return false;
    }

    public Object getDataObject(int row, int col)
    {
        return gridWeights[row][col][0].dataPoint.object;
    }
    
    public void addDataPoint(Object object, int row, int col)
    {
//		checkInitialize();
        if(nGridWeights[row][col] == -1) return;
        DataPoint point = new DataPoint(object, row, col);
        newDataPoints[row].add(point);
        DistanceSq weight = new DistanceSq(row, col, point, 0.0f);
        gridWeights[row][col] = new DistanceSq[]{weight};
        nGridWeights[row][col] = -1;
        setGridStatus(row, col, 0.0f);
        upToDate = false;
        nTotalPoints++;
    }

    public void deleteDataPoint(Object object, int row, int col)
    {
        nGridWeights[row][col] = 0;
        gridWeights[row][col] = null;
        isNeighbor[row][col] = false;
        gridChanged[row][col] = true;
        rowChanged[row] = true;
        if(colChanged != null) colChanged[col] = true;
        upToDate = false;
        nTotalPoints--;
    }

    protected void moveNewDataPoints()
    {
        for (int row = 0; row < nRows; row++)
            moveNewDataPoints(row);
    }

    protected void moveNewDataPoints(int row)
    {
        int nDataPoints = newDataPoints[row].size();
        for (int n = 0; n < nDataPoints; n++)
            dataPoints[row].add(newDataPoints[row].get(n));
        newDataPoints[row].clear();
    }

    /** if this point is interpolated from this object, set its grid status as changed */
    protected boolean updatePoint(Object object, int row, int col)
    {
        try
        {
            int nWeights = nGridWeights[row][col];
            if(nWeights <= 0) return false; // this is a data point
            DistanceSq[] pointWeights = gridWeights[row][col];
            for(int n = 0; n < nWeights; n++)
            {
                DistanceSq pointWeight = pointWeights[n];
                if(pointWeight.dataPoint.object == object)
                {
    //				if(mainDebug) System.out.println("Updated row " + row + " col " + col);
                    setGridStatus(row, col, pointWeight.distanceSq);
                    nUpdatedPoints++;
                    return true;
                }
            }
            return false;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "updatePoint", " Failed at row: " + row + " col: " + col, e);
            return false;
        }
    }

    class DataPoint
    {
        Object object;
        int row, col;

        DataPoint(Object object, int row, int col)
        {
            this.object = object;
            this.row = row;
            this.col = col;
        }
    }

    class DistanceSq
    {
        DataPoint dataPoint;
        double distanceSq;

        DistanceSq(int row, int col, DataPoint dataPoint, double distanceSq)
        {
            this.dataPoint = dataPoint;
            this.distanceSq = distanceSq;
        }

        void print()
        {
            System.out.print(" row " + dataPoint.row + " col " + dataPoint.col + " distSq " + distanceSq);
        }
    }

    public class Weights
    {
        public int nWeights;
        public Object[] dataObjects;
        public double[] weights;
        public double weightSum;

        Weights(int nWeights, DistanceSq[] distanceSqs)
        {
            try
            {
                if (nWeights == -1 || nWeights == 1)
                {
                    this.nWeights = 1;
                    weights = new double[1];
                    dataObjects = new Object[1];
                    dataObjects[0] = distanceSqs[0].dataPoint.object;
                    weights[0] = 1.0;
                    weightSum = 1.0;
                    return;
                }
                if (nWeights < nMaxPoints || !modify)
                {
                    this.nWeights = nWeights;
                    dataObjects = new Object[nWeights];
                    double[] distances = new double[nWeights];
                    for (int n = 0; n < nWeights; n++)
                    {
                        dataObjects[n] = distanceSqs[n].dataPoint.object;
                        distances[n] = distanceSqs[n].distanceSq;
//						distances[n] = Math.sqrt(distanceSqs[n].distanceSq);
                    }
                    weightSum = 0.0;
                    weights = new double[nWeights];
                    for (int n = 0; n < nWeights; n++)
                    {
                        weights[n] = (distances[n] == 0.0) ? 1.0 : 1.0 / distances[n];
                        weightSum += weights[n];
                    }
                    if (weightSum == 0.0)
                    {
//					System.out.println("StsSpiralRadialInterpolation.Weights failed. weightSum = 0.0");
                    }

                }
                else
                {
                    nWeights--; // last object will be weighted zero, so we don't need it; just it's distance
                    this.nWeights = nWeights;
                    dataObjects = new Object[nWeights];
                    double[] distances = new double[nWeights];
                    for (int n = 0; n < nWeights; n++)
                    {
                        dataObjects[n] = distanceSqs[n].dataPoint.object;
                        distances[n] = Math.sqrt(distanceSqs[n].distanceSq);
                    }
                    double maxDistance = Math.sqrt(distanceSqs[nWeights].distanceSq);
                    weightSum = 0.0;
                    weights = new double[nWeights];
                    for (int n = 0; n < nWeights; n++)
                    {
                        double dwt = (distances[n] - maxDistance);
                        weights[n] = dwt * dwt / (distances[n] * maxDistance);
                        weightSum += weights[n];
                    }
                    if (weightSum == 0.0)
                    {
//					System.out.println("StsSpiralRadialInterpolation.Weights failed. weightSum = 0.0");
                    }
                }
            }
            catch (Exception e)
            {
                StsException.outputException("StsSpiralRadialInterpolation.Weights() failed.", e, StsException.WARNING);
            }
        }
    }
}