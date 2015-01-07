package com.Sts.Utilities.Interpolation;


import com.Sts.DBTypes.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

import java.util.*;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */

/** Does 1/R interpolation on a regular grid.  If the number of neighbor points is greater than the number of max points desired,
 *  a modified interpolation is used where weights scale to zero at the furthest point to prevent jumps in the interpolation as
 *  you move between regions.  
 */
public class StsSpiralRadialInterpolation extends StsRadialInterpolation
{
    float dRow;
	float dCol;
	float dRowSq;
	float dColSq;
    float[] rowDSqs;
	float[] colDSqs;
    int maxNSpirals;

	public StsSpiralRadialInterpolation(int nRows, int nCols, float dRow, float dCol, int nMinPoints, int nMaxPoints, float neighborRadius)
	{
		this.nRows = nRows;
		this.nCols = nCols;
		this.dRow = dRow;
		this.dCol = dCol;
		this.dRowSq = dRow*dRow;
		this.dColSq = dCol*dCol;
		this.nMinPoints = nMinPoints;
		this.nMaxPoints = nMaxPoints;
		this.neighborRadius = neighborRadius;
		this.maxDistSq = neighborRadius*neighborRadius;
		this.maxNSpirals = Math.max(nRows, nCols);

//		classInitialize();
	}

    public StsSpiralRadialInterpolation(StsPreStackLineSet3d seismicVolume, int nMinPoints, int nMaxPoints, float neighborRadius)
	{
        this.lineSet = seismicVolume;
        this.nRows = seismicVolume.nRows;
		this.nCols = seismicVolume.nCols;
		this.dRow = seismicVolume.yInc;
		this.dCol = seismicVolume.xInc;
		this.dRowSq = dRow*dRow;
		this.dColSq = dCol*dCol;
		this.nMinPoints = nMinPoints;
		this.nMaxPoints = nMaxPoints;
		this.neighborRadius = neighborRadius;
		this.maxDistSq = neighborRadius*neighborRadius;
		this.maxNSpirals = Math.max(nRows, nCols);
        initialize();
    }

    public StsSpiralRadialInterpolation(StsRotatedGridBoundingBox boundingBox, int nMinPoints, int nMaxPoints)
	{
		this.nRows = boundingBox.nRows;
		this.nCols = boundingBox.nCols;
		this.dRow = boundingBox.yInc;
		this.dCol = boundingBox.xInc;
		this.dRowSq = dRow*dRow;
		this.dColSq = dCol*dCol;
		this.nMinPoints = nMinPoints;
		this.nMaxPoints = nMaxPoints;
		this.neighborRadius = StsParameters.largeFloat;
		this.maxDistSq = neighborRadius*neighborRadius;
		this.maxNSpirals = Math.max(nRows, nCols);
        initialize();
    }

	private void checkInitialize()
	{
		if(gridWeights != null) return;
		initialize();
	}

	public void initialize()
	{
		gridWeights = new DistanceSq[nRows][nCols][];
		nGridWeights = new int[nRows][nCols];
		rowDSqs = new float[nRows+2];
		colDSqs = new float[nCols+2];
//		nOldGridWeights = nGridWeights;
//		nGridWeights = new int[nRows][nCols];
		gridChanged = new boolean[nRows][nCols];
		isNeighbor = new boolean[nRows][nCols];
		rowChanged = new boolean[nRows];
		colChanged = new boolean[nCols];
        dataPoints = new ArrayList[nRows];
        newDataPoints = new ArrayList[nRows];
        nTotalPoints = 0;
        for(int row = 0; row < nRows; row++)
        {
            dataPoints[row] = new ArrayList();
            newDataPoints[row] = new ArrayList();
        }
        upToDate = false;
    }

	public void run()
	{
        if(nTotalPoints == 0) return;
        int nSpiral = 0;
		boolean runSpiral = true;
		while(runSpiral && nSpiral < maxNSpirals)
		{
			runSpiral = false;
			nSpiral++;
			currentMinNWeights = StsParameters.largeInt;
			if(!computeDistanceSquares(nSpiral)) break; // break if spiral is beyond maxRadius
            for(int row = 0; row < nRows; row++)
            {
                Object[] points = newDataPoints[row].toArray();
		        int nPoints = points.length;
                for(int n = 0; n < nPoints; n++)
                    if(runPointSpiral((DataPoint)points[n], nSpiral)) runSpiral = true;
            }
            if(debug) System.out.println("Spiral " + nSpiral + " currentMinNWeights " + currentMinNWeights);
		}
		if(debug) System.out.println("Last spiral " + nSpiral);
		upToDate = true;
		moveNewDataPoints();
	}

    private boolean computeDistanceSquares(int nSpiral)
	{
		rowDSqs[0] = dColSq*nSpiral*nSpiral;
		colDSqs[0] = dRowSq*nSpiral*nSpiral;
//		if(rowDSqs[0] > maxDistSq && colDSqs[0] > maxDistSq) return false;

		int nn = Math.min(nRows, nSpiral) + 1;
		for(int n = 1; n < nn; n++)
			rowDSqs[n] = rowDSqs[0] + dRowSq * n * n;
		nn = Math.min(nCols, nSpiral) + 1;
		for(int n = 1; n < nn; n++)
			colDSqs[n] = colDSqs[0] + dColSq * n * n;
		return true;
	}

	private boolean runPointSpiral(DataPoint dataPoint, int nSpiral)
	{
		int rowMin = dataPoint.row - nSpiral;
		int rowMax = dataPoint.row + nSpiral;
		int colMin = dataPoint.col - nSpiral;
		int colMax = dataPoint.col + nSpiral;

		int rowCenter = dataPoint.row;
		int colCenter = dataPoint.col;

		boolean terminate = true;
		boolean spiralOK = false;
		boolean doBottom = rowMin >= 0;
		if(doBottom)
		{
			int lastCol = Math.min(colMax-1, nCols-1);
			for(int col = colCenter, n = 0; col <= lastCol; col++, n++)
				if(addWeight(dataPoint, rowMin, col, rowDSqs[n])) spiralOK = true;
			lastCol = Math.max(0, colMin);
			for(int col = colCenter-1, n = 1; col >= lastCol; col--, n++)
				if(addWeight(dataPoint, rowMin, col, rowDSqs[n])) spiralOK = true;
		}

		boolean doTop = rowMax < nRows;
		if(doTop)
		{
			int lastCol = Math.min(colMax, nCols-1);
			for(int col = colCenter+1, n = 1; col <= lastCol; col++, n++)
				if(addWeight(dataPoint, rowMax, col, rowDSqs[n])) spiralOK = true;
			lastCol = Math.max(0, colMin+1);
			for(int col = colCenter, n = 0; col >= lastCol; col--, n++)
				if(addWeight(dataPoint, rowMax, col, rowDSqs[n])) spiralOK = true;
		}

		boolean doRight = colMax < nCols;
		if(doRight)
		{
			int lastRow = Math.min(rowMax-1, nRows-1);
			for(int row = rowCenter, n = 0; row <= lastRow; row++, n++)
				if(addWeight(dataPoint, row, colMax, rowDSqs[n])) spiralOK = true;
			lastRow = Math.max(0, rowMin);
			for(int row = rowCenter-1, n = 1; row >= lastRow; row--, n++)
				if(addWeight(dataPoint, row, colMax, rowDSqs[n])) spiralOK = true;
		}
		boolean doLeft = colMin >= 0;
		if(doLeft)
		{
			int lastRow = Math.min(rowMax, nRows-1);
			for(int row = rowCenter+1, n = 1; row <= lastRow; row++, n++)
				if(addWeight(dataPoint, row, colMin, rowDSqs[n])) spiralOK = true;
			lastRow = Math.max(0, rowMin+1);
			for(int row = rowCenter, n = 0; row >= lastRow; row--, n++)
				if(addWeight(dataPoint, row, colMin, rowDSqs[n])) spiralOK = true;
		}
		return spiralOK;
	}

    /** spirals out from existing point, flagging each point interpolated from this one as having been changed */
    public void runUpdate(Object object, int row, int col)
	{
        int nSpiral = 0;
		boolean runSpiral = true;
		this.nUpdatedPoints = 0;
		while(runSpiral && nSpiral < maxNSpirals)
		{
			nSpiral++;
			runSpiral = false;
			if(update(object, row, col, nSpiral)) runSpiral = true;
		}
        if(debug) System.out.println("Last spiral " + nSpiral + " nUpdatedPoints " + nUpdatedPoints);
    }

    /** runs one spiral about this existing point, flagging each point interpolated from this one as having been changed */
    private boolean update(Object object, int rowCenter, int colCenter, int nSpiral)
	{
		int rowMin = rowCenter - nSpiral;
		int rowMax = rowCenter + nSpiral;
		int colMin = colCenter - nSpiral;
		int colMax = colCenter + nSpiral;

		boolean terminate = true;
		boolean spiralOK = false;
		boolean doBottom = rowMin >= 0;
		if(doBottom)
		{
			int lastCol = Math.min(colMax-1, nCols-1);
			for(int col = colCenter, n = 0; col <= lastCol; col++, n++)
				if(updatePoint(object, rowMin, col)) spiralOK = true;
			lastCol = Math.max(0, colMin);
			for(int col = colCenter-1, n = 1; col >= lastCol; col--, n++)
				if(updatePoint(object, rowMin, col)) spiralOK = true;
		}

		boolean doTop = rowMax < nRows;
		if(doTop)
		{
			int lastCol = Math.min(colMax, nCols-1);
			for(int col = colCenter+1, n = 1; col <= lastCol; col++, n++)
				if(updatePoint(object, rowMax, col)) spiralOK = true;
			lastCol = Math.max(0, colMin+1);
			for(int col = colCenter, n = 0; col >= lastCol; col--, n++)
				if(updatePoint(object, rowMax, col)) spiralOK = true;
		}

		boolean doRight = colMax < nCols;
		if(doRight)
		{
			int lastRow = Math.min(rowMax-1, nRows-1);
			for(int row = rowCenter, n = 0; row <= lastRow; row++, n++)
				if(updatePoint(object, row, colMax)) spiralOK = true;
			lastRow = Math.max(0, rowMin);
			for(int row = rowCenter-1, n = 1; row >= lastRow; row--, n++)
				if(updatePoint(object, row, colMax)) spiralOK = true;
		}
		boolean doLeft = colMin >= 0;
		if(doLeft)
		{
			int lastRow = Math.min(rowMax, nRows-1);
			for(int row = rowCenter+1, n = 1; row <= lastRow; row++, n++)
				if(updatePoint(object, row, colMin)) spiralOK = true;
			lastRow = Math.max(0, rowMin+1);
			for(int row = rowCenter, n = 0; row >= lastRow; row--, n++)
				if(updatePoint(object, row, colMin)) spiralOK = true;
		}
		return spiralOK;
	}

    public boolean isPlaneOK(int dir, int nPlane)
	{
		switch(dir)
		{
			case StsParameters.XDIR:
				return !colChanged[nPlane];
			case StsParameters.YDIR:
				return !rowChanged[nPlane];
			default:
				return true;
		}
	}

	public void print()
	{
		for(int row = 0; row < nRows; row++)
			for(int col = 0; col < nCols; col++)
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
		StsSpiralRadialInterpolation interpolation = new StsSpiralRadialInterpolation(10, 10, 1.0f, 1.0f, 1, 3, 10.0f);
        interpolation.checkInitialize();
        interpolation.addDataPoint(null, 5, 5);
		interpolation.addDataPoint(null, 4, 4);
		interpolation.addDataPoint(null, 8, 8);
		interpolation.addDataPoint(null, 8, 2);
		interpolation.addDataPoint(null, 2, 8);
		interpolation.run();
		interpolation.print();
		System.out.println("\nInserting new point.\n");
	    interpolation.addDataPoint(null, 4, 5);
		interpolation.run();
		interpolation.print();
	}

}
