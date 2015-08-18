package com.Sts.Utilities;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version beta 1.0
 */

import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;

import java.util.*;

public class StsDistanceTransform
{
    public float[][] distances;
    int nRows, nCols;
    float rowInc, colInc;
    boolean converged;

    transient float debugMaxSweepDistance;
    transient float debugMaxChange;

    static final boolean debugDistTrans = false;
    static final int debugDistTransRow = 258;
    static final int debugDistTransCol = 437;

    static final float largeFloat = StsParameters.largeFloat;

    public StsDistanceTransform(StsRotatedGridBoundingBox boundingBox)
    {
        this.nRows = boundingBox.getNRows();
        this.nCols = boundingBox.getNCols();
        this.rowInc = boundingBox.yInc;
        this.colInc = boundingBox.xInc;
        StsDistanceTransformPoint.setNRowsCols(nRows, nCols);
        initializeDistances();
    }

    public void distanceTransform(StsProgressPanel panel)
    {
        int row, col;
        float distance;

        try
        {
            if (distances == null) return;
            StsDistanceTransformPoint[] points = new StsDistanceTransformPoint[5];
            float diagonalDistance = (float) Math.sqrt(rowInc * rowInc + colInc * colInc);
            points[0] = new StsDistanceTransformPoint(0f);
            points[1] = new StsDistanceTransformPoint(diagonalDistance);
            points[2] = new StsDistanceTransformPoint(colInc);
            points[3] = new StsDistanceTransformPoint(diagonalDistance);
            points[4] = new StsDistanceTransformPoint(rowInc);

            converged = false;
            int iter = 0;
            int maxIter = 5;
            int nPoints = nRows * nCols;
            // max number of progress bar steps is maxIter times number of sweeps per iteration (4)
            if (panel != null)
            {
                panel.progressBar.setMaximum(maxIter*4);
                panel.appendLine("Interpolating distance using distance transform interpolation....");
            }
            long count = 0;
            while (!converged)
            {
                converged = true;
                // lower left to upper right sweep
                // for a point p0 computed weighted average value and nearest distance
                // from lower left 4 points and center point:
                //
                // p4 p0
                // p1 p2 p3
                //
                if (debugDistTrans) initializeDebugSweep(iter, "Sweep LL to UR");
                for (row = 0; row < nRows; row++)
                {
                    for (col = 0; col < nCols; col++)
                    {
                        distance = points[0].initialize(row, col, distances);
                        if (distance != 0.0f)
                        {
                            points[1].initialize(row - 1, col - 1, distances);
                            points[2].initialize(row - 1, col, distances);
                            points[3].initialize(row - 1, col + 1, distances);
                            points[4].initialize(row, col - 1, distances);
                            distanceTransform(row, col, points);
                        }
                    }
                }
                if(debugDistTrans) debugPrintDistance();
                if (panel != null)
                {
    //                panel.appendLine("Distance transform interpolation sweep Sweep LL to UR");
                    count += nPoints;
                    panel.progressBar.incrementProgress();
                }
                // upper right to lower left sweep
                // for a point p0 computed weighted average value and nearest distance
                // from upper left 4 points:
                //
                // p1 p2 p3
                //    p0 p4
                //

                if (debugDistTrans) initializeDebugSweep(iter, "Sweep UR to LL");
                for (row = nRows - 1; row >= 0; row--)
                {
                    for (col = nCols - 1; col >= 0; col--)
                    {
                        distance = points[0].initialize(row, col, distances);
                        if (distance != 0.0f)
                        {
                            points[1].initialize(row + 1, col - 1, distances);
                            points[2].initialize(row + 1, col, distances);
                            points[3].initialize(row + 1, col + 1, distances);
                            points[4].initialize(row, col + 1, distances);
                            distanceTransform(row, col, points);
                        }
                    }
                }
                if(debugDistTrans) debugPrintDistance();
                if (panel != null)
                {
    //                panel.appendLine("Distance transform interpolation sweep Sweep UR to LL");
                    count += nPoints;
                    panel.progressBar.incrementProgress();
                }
                // lower right to upper left sweep
                // for a point p0 computed weighted average value and nearest distance
                // from lower left 4 points and center point:
                //
                //    p0 p4
                // p1 p2 p3

                if (debugDistTrans) initializeDebugSweep(iter, "Sweep LR to UL");
                for (row = 0; row < nRows; row++)
                {
                    for (col = nCols - 1; col >= 0; col--)
                    {
                        distance = points[0].initialize(row, col, distances);
                        if (distance != 0.0f)
                        {
                            points[1].initialize(row - 1, col - 1, distances);
                            points[2].initialize(row - 1, col, distances);
                            points[3].initialize(row - 1, col + 1, distances);
                            points[4].initialize(row, col + 1, distances);
                            distanceTransform(row, col, points);
                        }
                    }
                }
                if(debugDistTrans) debugPrintDistance();
                if (panel != null)
                {
    //                panel.appendLine("Distance transform interpolation sweep Sweep LR to UL");
                    count += nPoints;
                    panel.progressBar.incrementProgress();
                }
                // upper left to lower right sweep
                // for a point p0 computed weighted average value and nearest distance
                // from upper left 4 points:
                //
                // p1 p2 p3
                // p4 p0

                if (debugDistTrans) initializeDebugSweep(iter, "Sweep UL to LR");
                for (row = nRows - 1; row >= 0; row--)
                {
                    for (col = 0; col < nCols; col++)
                    {
                        distance = points[0].initialize(row, col, distances);
                        if (distance != 0.0f)
                        {
                            points[1].initialize(row + 1, col - 1, distances);
                            points[2].initialize(row + 1, col, distances);
                            points[3].initialize(row + 1, col + 1, distances);
                            points[4].initialize(row, col - 1, distances);
                            distanceTransform(row, col, points);
                        }
                    }
                }
                if(debugDistTrans) debugPrintDistance();
                if (panel != null)
                {
    //                panel.appendLine("Distance transform interpolation sweep Sweep UL to LR");
                    count += nPoints;
                    panel.progressBar.incrementProgress();
                }
                iter++;
                if (!converged && iter >= maxIter)
                {
                    float largestDistance = 0.0f;
                    int largestDistanceRow = -1;
                    int largestDistanceCol = -1;
                    col = 0;
                    for (row = 0; row < nRows; row++)
                    {
                        for (col = 0; col < nCols; col++)
                        {
                            distance = distances[row][col];
                            if (distance == largeFloat)
                            {
                                StsException.systemError(this, "distanceTransform",  "failed to converge after " + maxIter + " iterations.\n" +
                                    "    First null value at row: " + row + " col: " + col);
                                return;
                            }
                            else if(distance > largestDistance)
                            {
                                largestDistance = distance;
                                largestDistanceRow = row;
                                largestDistanceCol = col;
                            }
                        }
                    }
                    StsException.systemError(this, "distanceTransform",  "failed to converge after " + maxIter + " iterations.\n" +
                                    "    Largest distance value is " + largestDistance + " at row: " + largestDistanceRow + " col: " + largestDistanceCol);
                    return;
                }
                else if (converged)
                    StsMessageFiles.logMessage("Distance transform converged in " + iter + " iterations.");
            }
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "distanceTransform", e); 
        }
    }

    private void initializeDebugSweep(int iter, String sweepName)
    {
        debugMaxSweepDistance = 0.0f;
        debugMaxChange = 0.0f;
        System.out.print("Iteration " + iter + " " + sweepName);
    }

    private void debugPrintDistance()
    {
        System.out.println(" maxDistance: " + debugMaxSweepDistance + " maxChange: " + debugMaxChange);
    }

    private void initializeDistances()
    {
        distances = new float[nRows][nCols];
        for (int row = 0; row < nRows; row++)
            Arrays.fill(distances[row], largeFloat);
    }

    public void setDistance(int row, int col, float distance)
    {
        if(distance < distances[row][col])
            distances[row][col] = distance;
    }

    public void distanceTransform(int row, int col, StsDistanceTransformPoint[] points)
    {
        float currentDistance = points[0].distance;
        float nearestDistance = largeFloat;
        boolean rowColDebug = false;
        if(debugDistTrans)
        {
            rowColDebug = debugDistanceTransform(row, col);
            if(rowColDebug)
                System.out.println("    row: " + row + " col: " + col + " currentDistance: " + currentDistance);
        }
        for (int n = 1; n < 5; n++)
        {
            float distance = points[n].getDistance();
            if (distance < nearestDistance)
            {
                nearestDistance = distance;
            }
        }
        if (nearestDistance < currentDistance)
        {
            converged = false;
            distances[row][col] = nearestDistance;
            if(debugDistTrans)
            {
                debugMaxSweepDistance = Math.max(nearestDistance, debugMaxSweepDistance);
                if(currentDistance != largeFloat)
                    debugMaxChange = Math.max(debugMaxChange, currentDistance - nearestDistance);
                if(rowColDebug)
                {
                    System.out.println("    new distance set: " + nearestDistance);
                    return;
                }
            }
        }
        else if(currentDistance == largeFloat)
        {
            converged = false;
        }
        if(rowColDebug)
        {
            System.out.println();
        }
    }

    boolean debugDistanceTransform(int row, int col)
    {
        return row == debugDistTransRow && col == debugDistTransCol;
    }
}