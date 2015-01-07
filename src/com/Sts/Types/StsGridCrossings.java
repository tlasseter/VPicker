
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.Interfaces.*;
import com.Sts.Utilities.*;

import java.util.*;

/** Generates a sequence of StsGridCrossingPoints which have
 *  index coordinates, rowOrCol, connect direction, and interpolation factor
 *  for any grid crossings between two points
 */
public class StsGridCrossings
{
    public ArrayList gridPoints = null;  /** List of grid crossing points */

    protected static final float largeFloat = StsParameters.largeFloat;

    protected static final int MINUS = StsParameters.MINUS;
    protected static final int PLUS = StsParameters.PLUS;
    protected static final int NONE = StsParameters.NONE;

    protected static final int ROW = StsParameters.ROW;
    protected static final int COL = StsParameters.COL;
    protected static final int ROWCOL = StsParameters.ROWCOL;

    static float gridLength = 0.0f;

    public StsGridCrossings(StsPoint point0, StsPoint point1, StsXYSurfaceGridable grid)
    {
		this(new StsGridCrossingPoint(grid, point0), new StsGridCrossingPoint(grid, point1), true);
    }

    /** use this constructor where we are interested in grid intersection of XY line between points */
    public StsGridCrossings(StsPoint point0, StsPoint point1, StsXYGridable grid)
    {
		this(new StsGridCrossingPoint(grid, point0), new StsGridCrossingPoint(grid, point1), true);
    }

    static public ArrayList<StsGridCrossingPoint> computeGridCrossingPoints(StsPoint[] points, StsXYGridable grid)
    {
        int nPoints = points.length;
        StsGridCrossingPoint[] gridPoints = new StsGridCrossingPoint[nPoints];
        for(int n = 0; n < nPoints; n++)
            gridPoints[n] = new StsGridCrossingPoint(grid, points[n]);
        return computeGridCrossingPoints(gridPoints, grid);
    }

    static public ArrayList<StsGridCrossingPoint> computeGridCrossingPoints(StsGridCrossingPoint[] gridPoints, StsXYGridable grid)
    {
        ArrayList<StsGridCrossingPoint> gridCrossingsList = new ArrayList<StsGridCrossingPoint>();
        int nPoints = gridPoints.length;
        for(int n = 0; n < nPoints-1; n++)
        {
            StsGridCrossings gridCrossings = new StsGridCrossings(gridPoints[n], gridPoints[n+1], false);
            gridCrossingsList.addAll(gridCrossings.gridPoints);
        }
        gridCrossingsList.add(gridPoints[nPoints-1]);
        return gridCrossingsList;
    }

    public StsGridCrossings(StsGridCrossingPoint gridPoint0, StsGridCrossingPoint gridPoint1, boolean includeLastPoint)
    {
        int iStart = 0, iEnd = 0;               /** grid I lines in between the two points */
        int jStart = 0, jEnd = 0;               /** grid J lines in between the two points */
        float dJdI = 0.0f, dIdJ = 0.0f;         /** slope values used in computing intermediate values */
        float gridRowF, gridColF;               /** float IJ values at the intermediate point */
        int iConnect = MINUS;                   /** connect directions along J grid line */
        int  jConnect= MINUS;                   /** connect directions along J grid line */
        float iff = 0.0f, jff = 0.0f;           /** interpolation factors along grid lines */
        StsGridCrossingPoint iGridPoint = null;  /** gridPoint along I line */
        StsGridCrossingPoint jGridPoint = null;  /** gridPoint along J line */
        StsGridCrossingPoint gridCrossingPoint;          /** scratch gridPoint      */
        StsPoint point;
        boolean iOK, jOK;
        int iInc = 1, jInc = 1;
        int i, j;
        float dI, dJ;

        gridPoints = new ArrayList(10); // List of new points added

        float gridRowF0 = gridPoint0.iF;
        float gridColF0 = gridPoint0.jF;
        float gridRowF1 = gridPoint1.iF;
        float gridColF1 = gridPoint1.jF;

        // if points are in same cell, return the two points
        if((int)gridRowF1 == (int)gridRowF0 && (int)gridColF1 == (int)gridColF0)
        {
            gridPoints.add(gridPoint0);
            if(includeLastPoint) gridPoints.add(gridPoint1);
            return;
        }

        gridPoints.add(gridPoint0);

        if(gridRowF0 == gridRowF1)
            iOK = false;
        else
        {
            dJdI = (gridColF1 - gridColF0)/(gridRowF1 - gridRowF0);

            if(gridRowF1 > gridRowF0)
            {
                iStart = StsMath.ceiling(gridRowF0);
                iEnd =   StsMath.floor(gridRowF1);

                if(iStart == gridRowF0) iStart++;
                if(iEnd == gridRowF1) iEnd--;
                iInc = 1;
                iConnect = PLUS; /** We are going in +I direction: +J is to rite	*/
                iOK = iEnd >= iStart;
            }
            else // gridRowF1 < gridRowF0
            {
                iStart = StsMath.floor(gridRowF0);
                iEnd = StsMath.ceiling(gridRowF1);

                if(iStart == gridRowF0) iStart--;
                if(iEnd == gridRowF1) iEnd++;
                iInc = -1;
                iConnect = MINUS; /** We are going in -I direction: -J is to rite	*/
                iOK = iEnd <= iStart;
            }
        }

        if(gridColF0 == gridColF1)
            jOK = false;
        else
        {
            dIdJ = (gridRowF1 - gridRowF0)/(gridColF1 - gridColF0);

            if(gridColF1 > gridColF0)
            {
                jStart = StsMath.ceiling(gridColF0);
                jEnd = StsMath.floor(gridColF1);

                if(jStart == gridColF0) jStart++;
                if(jEnd == gridColF1) jEnd--;
                jInc = 1;
                jConnect = MINUS;  /** We are going in +J direction: -I is to rite	*/
                jOK = jEnd >= jStart;
            }
            else // gridColF1 < gridColF0
            {
                jStart = StsMath.floor(gridColF0);
                jEnd = StsMath.ceiling(gridColF1);

                if(jStart == gridColF0) jStart--;
                if(jEnd == gridColF1) jEnd++;
                jInc = -1;
                jConnect = PLUS;  /** We are going in -J direction: +I is to rite	*/
                jOK = jEnd <= jStart;
            }
        }

        gridPoint0.setRowOrCol();
        gridPoint0.setConnects(iConnect, jConnect);
        gridPoint1.setRowOrCol();
        gridPoint1.setConnects(iConnect, jConnect);

        i = iStart;
        j = jStart;

        /** compute location of next point on I grid lines */
        if(iOK)
        {
            gridRowF = (float)i;
            gridColF = gridColF0 + (i - gridRowF0)*dJdI;
            iff = (gridRowF - gridRowF0)/(gridRowF1 - gridRowF0);
            iGridPoint = new StsGridCrossingPoint(gridRowF, ROW, gridPoint0, gridPoint1, iff, iConnect);
        }
        else
            iff = largeFloat;

        /** compute location of next point on J grid lines */
        if(jOK)
        {
            gridColF = (float)j;
            jff = (gridColF - gridColF0)/(gridColF1 - gridColF0);
            jGridPoint = new StsGridCrossingPoint(gridColF, COL, gridPoint0, gridPoint1, jff, jConnect);
        }
        else
            jff = largeFloat;


        while(iOK || jOK)
        {
            if(iff <= jff)
            {
                gridPoints.add(iGridPoint);
                i = i + iInc;
                if(iInc > 0 && i <= iEnd || iInc < 0 && i >= iEnd)
                {
                    gridRowF = (float)i;
                    iff = (gridRowF - gridRowF0)/(gridRowF1 - gridRowF0);
                    iGridPoint = new StsGridCrossingPoint(gridRowF, ROW, gridPoint0, gridPoint1, iff, iConnect);
                }
                else
                {
                    iff = largeFloat;
                    iOK = false;
                }
            }
            else if(jff <= iff)
            {
                gridPoints.add(jGridPoint);
                j = j + jInc;
                if(jInc > 0 && j <= jEnd || jInc < 0 && j >= jEnd)
                {
                    gridColF = (float)j;
                    jff = (gridColF - gridColF0)/(gridColF1 - gridColF0);
                    jGridPoint = new StsGridCrossingPoint(gridColF, COL, gridPoint0, gridPoint1, jff, jConnect);
                }
                else
                {
                    jff = largeFloat;
                    jOK = false;
                }
            }
        }
        if(includeLastPoint) gridPoints.add(gridPoint1);
    }

	public StsGridPoint getGridIntersection(StsXYSurfaceGridable grid)
	{
		StsGridCrossingPoint gridPoint0, gridPoint1;

		if(gridPoints == null) return null;
		int nGridPoints = gridPoints.size();
		gridPoint1 = (StsGridCrossingPoint)gridPoints.get(0);
		for(int n = 1; n < nGridPoints; n++)
		{
		    gridPoint0 = gridPoint1;
			gridPoint1 = (StsGridCrossingPoint)gridPoints.get(n);
			StsGridCrossingPoint gridPoint = getGridIntersection(gridPoint0, gridPoint1, grid);
			if(gridPoint != null)
            {
                return new StsGridPoint(gridPoint);
            }
		}
		return null;
	}

	private StsGridCrossingPoint getGridIntersection(StsGridCrossingPoint gridPoint0, StsGridCrossingPoint gridPoint1, StsXYSurfaceGridable grid)
	{
		float z0 = gridPoint0.point.getZorT();
		float z1 = gridPoint1.point.getZorT();
		float surfZ0 = getGridZ(gridPoint0, grid);
        if(surfZ0 == StsParameters.nullValue) return null;
        float surfZ1 = getGridZ(gridPoint1, grid);
        if(surfZ1 == StsParameters.nullValue) return null;
        surfZ0 = StsMath.checkRoundOffInteger(surfZ0);
        surfZ1 = StsMath.checkRoundOffInteger(surfZ1);
        float dZ0 = z0 - surfZ0;
		if(dZ0 == 0.0f) return gridPoint0;
		float dZ1 = z1 - surfZ1;
		if(dZ1 == 0.0f) return gridPoint1;

        if (dZ0 <= 0.0f && dZ1 >= 0.0f || dZ0 >= 0.0f && dZ1 <= 0.0f )
        {
            float f = Math.abs(dZ0/(dZ1 - dZ0));
			return new StsGridCrossingPoint(gridPoint0, gridPoint1, f);
        }
		return null;
	}

	private float getGridZ(StsGridCrossingPoint gridPoint, StsXYSurfaceGridable grid)
	{
		StsGridPoint stsGridPoint = new StsGridPoint(gridPoint.iF, gridPoint.jF, grid);
        return grid.interpolateBilinearZ(stsGridPoint, true, true);
//		return grid.interpolateBilinearZ(stsGridPoint, true, false);
	}

    static public StsGridPoint[] getStsGridPoints(ArrayList gridCrossingPointsList)
    {
        int nPoints = gridCrossingPointsList.size();
        StsGridPoint[] gridPoints = new StsGridPoint[nPoints];
        for(int n = 0; n < nPoints; n++)
        {
            StsGridCrossingPoint gridPoint = (StsGridCrossingPoint) gridCrossingPointsList.get(n);
            gridPoints[n] = new StsGridPoint(gridPoint.point, gridPoint.iF, gridPoint.jF);
        }
        return gridPoints;
    }

    static public StsGridPoint[] getStsGridPoints(StsPoint[] points, StsXYGridable grid)
    {
        int nPoints = points.length;
        StsGridPoint[] gridPoints = new StsGridPoint[nPoints];
        for(int n = 0; n < nPoints; n++)
        {
            StsGridCrossingPoint gridPoint = new StsGridCrossingPoint(grid, points[n]);
            gridPoints[n] = new StsGridPoint(gridPoint.point, gridPoint.iF, gridPoint.jF);
        }
        return gridPoints;
    }
}
