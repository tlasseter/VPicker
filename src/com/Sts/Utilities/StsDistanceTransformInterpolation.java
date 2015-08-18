package com.Sts.Utilities;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version beta 1.0
 */

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;

public class StsDistanceTransformInterpolation
 {
    StsDistanceTransformFace distanceTransformObject;
    float[][] distances;
    int nRows, nCols;
    boolean debug = false;
    float maxInterpolationDistance;
    static final float maxInterpolationGridSpacing = 4.0f;

    static final float largeFloat = StsParameters.largeFloat;

    public StsDistanceTransformInterpolation(StsDistanceTransformFace distanceTransformObject, boolean debug)
    {
        this.distanceTransformObject = distanceTransformObject;
        this.nRows = distanceTransformObject.getNRows();
        this.nCols = distanceTransformObject.getNCols();
        StsDistanceTransformPoint.setNRowsCols(nRows, nCols);
        this.debug = debug;
    }

	public void interpolateDistanceTransform(StsProgressPanel panel)
	{
        int row, col;
        float distance, newDistance;

        float[][] distances = distanceTransformObject.initializeDistances();
        if(distances == null) return;
        StsDistanceTransformPoint[] points = new StsDistanceTransformPoint[5];
        float[] distanceParameters = distanceTransformObject.getDistanceParameters();
        float rowDistance = distanceParameters[0];
        float colDistance = distanceParameters[1];
        float diagonalDistance = (float)Math.sqrt(rowDistance*rowDistance + colDistance*colDistance);
        maxInterpolationDistance = maxInterpolationGridSpacing*diagonalDistance;
        points[0] = new StsDistanceTransformPoint(0f);
        points[1] = new StsDistanceTransformPoint(diagonalDistance);
        points[2] = new StsDistanceTransformPoint(colDistance);
        points[3] = new StsDistanceTransformPoint(diagonalDistance);
        points[4] = new StsDistanceTransformPoint(rowDistance);

        boolean converged = false;
        int iter = 0;
        int maxIter = 5;
        int intervalMin, intervalMax;
        int nPoints = nRows*nCols;
        // for the progress bar, we will assume that we will run about 2 iterations
        // or eight times thru row and col do loops
        if(panel != null)
        {
            panel.progressBar.setSubIntervalTotalCount(8*nRows*nCols);
            panel.appendLine("Interpolating null points using distance transform interpolation....");
        }
        long count = 0;
        while(!converged)
        {
            converged = true;
            // lower left to upper right sweep
            // for a point p0 computed weighted average value and nearest distance
            // from lower left 4 points and center point:
            //
            // p4 p0
            // p1 p2 p3
            //

            if(debug) System.out.println("Sweep LL to UR");
            for(row = 0; row < nRows; row++)
            {
                for(col = 0; col < nCols; col++)
                {
                    distance = points[0].initialize(row, col, distances);
                    if(distance != 0.0f)
                    {
                        points[1].initialize(row-1, col-1, distances);
                        points[2].initialize(row-1, col, distances);
                        points[3].initialize(row-1, col+1, distances);
                        points[4].initialize(row, col-1, distances);
                        newDistance = distanceTransformObject.doDistanceTransform(row, col, points);
                        if(newDistance != distance || newDistance == largeFloat) converged = false;
                        distances[row][col] = newDistance;
                    }
                }
            }
            if(panel != null)
            {
//                panel.appendLine("Distance transform interpolation sweep Sweep LL to UR");
                count += nPoints;
                panel.progressBar.setCount(count);
            }
            // upper right to lower left sweep
            // for a point p0 computed weighted average value and nearest distance
            // from upper left 4 points:
            //
            // p1 p2 p3
            //    p0 p4
            //

            if(debug) System.out.println("Sweep UR to LL");
            for(row = nRows-1; row >= 0; row--)
            {
                for(col = nCols-1; col >= 0; col--)
                {
                    distance = points[0].initialize(row, col, distances);
                    if(distance != 0.0f)
                    {
                        points[1].initialize(row+1, col-1, distances);
                        points[2].initialize(row+1, col, distances);
                        points[3].initialize(row+1, col+1, distances);
                        points[4].initialize(row, col+1, distances);
                        newDistance = distanceTransformObject.doDistanceTransform(row, col, points);
                        if(newDistance != distance || newDistance == largeFloat) converged = false;
                        distances[row][col] = newDistance;
                     }
                 }
            }
            if(panel != null)
            {
//                panel.appendLine("Distance transform interpolation sweep Sweep UR to LL");
                count += nPoints;
                panel.progressBar.setCount(count);
            }
            // lower right to upper left sweep
            // for a point p0 computed weighted average value and nearest distance
            // from lower left 4 points and center point:
            //
            //    p0 p4
            // p1 p2 p3

            if(debug) System.out.println("Sweep LR to UL");
            for(row = 0; row < nRows; row++)
            {
                for(col = nCols-1; col >= 0; col--)
                {
                    distance = points[0].initialize(row, col, distances);
                    if(distance != 0.0f)
                    {
                        points[1].initialize(row-1, col-1, distances);
                        points[2].initialize(row-1, col, distances);
                        points[3].initialize(row-1, col+1, distances);
                        points[4].initialize(row, col+1, distances);
                        newDistance = distanceTransformObject.doDistanceTransform(row, col, points);
                        if(newDistance != distance || newDistance == largeFloat) converged = false;
                        distances[row][col] = newDistance;
                     }
                }
            }
            if(panel != null)
            {
//                panel.appendLine("Distance transform interpolation sweep Sweep LR to UL");
                count += nPoints;
                panel.progressBar.setCount(count);
            }
            // upper left to lower right sweep
            // for a point p0 computed weighted average value and nearest distance
            // from upper left 4 points:
            //
            // p1 p2 p3
            // p4 p0

            if(debug) System.out.println("Sweep UL to LR");
            for(row = nRows-1; row >= 0; row--)
            {
                for(col = 0; col < nCols; col++)
                {
                    distance = points[0].initialize(row, col, distances);
                    if(distance != 0.0f)
                    {
                        points[1].initialize(row+1, col-1, distances);
                        points[2].initialize(row+1, col, distances);
                        points[3].initialize(row+1, col+1, distances);
                        points[4].initialize(row, col-1, distances);
                        newDistance = distanceTransformObject.doDistanceTransform(row, col, points);
                        if(newDistance != distance || newDistance == largeFloat) converged = false;
                        distances[row][col] = newDistance;
                     }
                }
            }
            if(panel != null)
            {
//                panel.appendLine("Distance transform interpolation sweep Sweep UL to LR");
                count += nPoints;
                panel.progressBar.setCount(count);
            }
            iter++;
            if(!converged && iter >= maxIter)
            {
                col = 0;
                for(row = 0; row < nRows; row++)
                    for(col = 0; col < nCols; col++)
                        if(distances[row][col] == largeFloat) break;
                StsException.systemError(this, "interpolateDistanceTransform",  "failed to converge after " + maxIter + " iterations.\n" +
                    "    First null value at row: " + row + " col: " + col);
                break;
            }
            else if(converged)
                StsMessageFiles.logMessage("Distance transform converged in " + iter + " iterations.");
        }
	}
}
