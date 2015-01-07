package com.Sts.Utilities.Interpolation;

import com.Sts.DBTypes.*;
import com.Sts.Utilities.*;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jan 29, 2007
 * Time: 3:34:15 PM
 * To change this template use File | Settings | File Templates.
 */

/** Does 1/R interpolation on a series of lines.  Each row is a line, the number of columns in each row is the number of points in the line.
 *  Each line is interpolated independently.
 */
public class StsLinesRadialInterpolation extends StsRadialInterpolation
{
    public float[][] pointX; // x coordinate of points on lines
    public float[][] pointY; // y coordinate of points on lines
    public StsPreStackLineSet2d seismicVolume = null;

    public StsLinesRadialInterpolation(StsPreStackLineSet2d seismicVolume, int nMinPoints, int nMaxPoints, float neighborRadius)
	{
        this.seismicVolume = seismicVolume;
        nRows = seismicVolume.nRows;
        this.nMinPoints = nMinPoints;
		this.nMaxPoints = nMaxPoints;
		this.neighborRadius = neighborRadius;
		this.maxDistSq = neighborRadius*neighborRadius;
        initialize();
    }

    public StsLinesRadialInterpolation(int nRows, int nMinPoints, int nMaxPoints, float neighborRadius)
	{
        this.nRows = nRows;
        this.nMinPoints = nMinPoints;
		this.nMaxPoints = nMaxPoints;
		this.neighborRadius = neighborRadius;
		this.maxDistSq = neighborRadius*neighborRadius;
        initialize();
    }

    public void initialize()
    {
        this.nColsInRow = new int[nRows];
        pointX = new float[nRows][];
        pointY = new float[nRows][];
        gridWeights = new DistanceSq[nRows][][];
        nGridWeights = new int[nRows][];
		isNeighbor = new boolean[nRows][];
		rowChanged = new boolean[nRows];
        gridChanged = new boolean[nRows][];
        dataPoints = new ArrayList[nRows];
        newDataPoints = new ArrayList[nRows];
        nTotalPoints = 0;
        for(int row = 0; row < nRows; row++)
        {
            dataPoints[row] = new ArrayList();
            newDataPoints[row] = new ArrayList();
        }

        if(seismicVolume == null) return;

        for(int row = 0; row < nRows; row++)
        {
            StsPreStackLine2d line = (StsPreStackLine2d)seismicVolume.lines[row];
            int nCols = line.nCols;
            this.nColsInRow[row] = nCols;
            this.pointX[row] = line.cdpX;
            this.pointY[row] = line.cdpY;
            gridWeights[row] = new DistanceSq[nCols][];
            nGridWeights[row] = new int[nCols];
            gridChanged[row] = new boolean[nCols];
            isNeighbor[row] = new boolean[nCols];
        }
        upToDate = false;
    }

    public void addLine(int row, int nCols, float[] pointX, float[] pointY)
    {
        this.nColsInRow[row] = nCols;
        this.pointX[row] = pointX;
        this.pointY[row] = pointY;
        gridWeights[row] = new DistanceSq[nCols][];
        nGridWeights[row] = new int[nCols];
        gridChanged[row] = new boolean[nCols];
        isNeighbor[row] = new boolean[nCols];
    }

    public void run()
	{
        if(nTotalPoints > 0)
        {
            for(int row = 0; row < nRows; row++)
                interpolateLine(row);
        }
        upToDate = true;
    }

    public void interpolateLine(int row)
    {
        int nStep = 0;
        int maxNSteps = nColsInRow[row];

        boolean runLine = true;
        while(runLine && nStep < maxNSteps)
        {
            runLine = false;
            nStep++;
            currentMinNWeights = StsParameters.largeInt;
            Object[] linePoints = newDataPoints[row].toArray();
            int nLinePoints = linePoints.length;
            for(int n = 0; n < nLinePoints; n++)
                if(runLineStep((DataPoint)linePoints[n], nStep)) runLine = true;
        }
 		if(debug) System.out.println("Last step " + nStep);
		moveNewDataPoints(row);
	}

    private boolean runLineStep(DataPoint dataPoint, int nStep)
	{
        double x, y, dx, dy, distSq;

        boolean terminate = true;
        boolean stepOK = false;
        int row = dataPoint.row;

        int colMin = dataPoint.col - nStep;
		int colMax = dataPoint.col + nStep;
        if(colMin < 0 && colMax >= nColsInRow[row]) return false;

        int rowCenter = dataPoint.row;
		int colCenter = dataPoint.col;

        double xCenter = pointX[rowCenter][colCenter];
        double yCenter = pointY[rowCenter][colCenter];

        if(colMin >= 0)
        {
            x = pointX[row][colMin];
            y = pointY[row][colMin];

            dx = x - xCenter;
            dy = y - yCenter;
            distSq = dx*dx + dy*dy;
            if(addWeight(dataPoint, row, colMin, distSq)) stepOK = true;
        }
        if(colMax < nColsInRow[row])
        {
            x = pointX[row][colMax];
            y = pointY[row][colMax];

            dx = x - xCenter;
            dy = y - yCenter;
            distSq = dx*dx + dy*dy;
            if(addWeight(dataPoint, row, colMax, distSq)) stepOK = true;
        }
		return stepOK;
	}

    public void runUpdate(Object object, int row, int col)
	{
        int nStep = 0;
        int maxNSteps = nColsInRow[row];
        boolean runLine = true;
        this.nUpdatedPoints = 0;

        while(runLine && nStep < maxNSteps)
        {
            runLine = false;
            nStep++;
            if(update(object, row, col, nStep)) runLine = true;
        }
 		dataPoints[row].add(object);
    }

	private boolean update(Object object, int rowCenter, int colCenter, int nStep)
	{
        boolean stepOK = false;

        int colMin = colCenter - nStep;
		int colMax = colCenter + nStep;
        if(colMin < 0 && colMax >= nColsInRow[rowCenter]) return false;

        if(colMin >= 0)
        {
            if(updatePoint(object, rowCenter, colMin)) stepOK = true;
        }
        if(colMax < nColsInRow[rowCenter])
        {
            if(updatePoint(object, rowCenter, colMax)) stepOK = true;
        }
		return stepOK;
    }

	public void print()
	{
		for(int row = 0; row < nRows; row++)
			for(int col = 0; col < nColsInRow[row]; col++)
			{
				System.out.print("row " + row + " col " + col + ":  ");
				int nGridWeights = this.nGridWeights[row][col];
				if(nGridWeights == -1) nGridWeights = 1;
				for(int n = 0; n < nGridWeights; n++)
				{
					DistanceSq weight = gridWeights[row][col][n];
					if(weight != null) weight.print();
				}
				System.out.println("");
			}
	}

	public static void main(String[] args)
	{
		StsLinesRadialInterpolation interpolation = new StsLinesRadialInterpolation(1, 1, 3, 10.0f);
        float[] pointsX = new float[10];
        float[] pointsY = new float[10];
        for(int n = 0; n < 10; n++)
        {
            pointsX[n] = n;
            pointsY[n] = 0;
        }
        interpolation.addLine(0, 10, pointsX, pointsY);

        interpolation.addDataPoint(null, 0, 3);
		interpolation.addDataPoint(null, 0, 4);
		interpolation.addDataPoint(null, 0, 7);
		interpolation.interpolateLine(0);
		interpolation.print();
//		System.out.println("\nInserting new point.\n");
//	    interpolation.addDataPoint(null, 4, 5);
//		interpolation.run();
//		interpolation.print();
	}
}
