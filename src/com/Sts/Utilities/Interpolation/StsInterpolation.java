package com.Sts.Utilities.Interpolation;

import com.Sts.MVC.*;
import com.Sts.Utilities.*;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jan 26, 2008
 * Time: 7:54:07 AM
 * To change this template use File | Settings | File Templates.
 */
public abstract class StsInterpolation
{
    static public StsInterpolation interpolator;
    /** horizon grid to be interpolated or for blockGrid interpolation, the grid to be interpolated from */
//    StsModelSurface surface;
    int iCenter; // row of point we want to interpolate
    int jCenter; // col of point we want to interpolate
    boolean useGradient; // use gradient in estimating point z
    public boolean debugInterpolate; // print mainDebug information
    boolean isConverged;        // indicates interpolation procedure converged
    int rowMin;
    int rowMax;
    int colMin;
    int colMax; // size of bounding box for interpolation
    boolean hasNulls = false;// parameters set once for all subsequent interpolations
    int nAvailablePoints = 0; // total number of good points existing on this grid
    int minNPoints = 10;
    int maxNPoints = 15;
    int maxNSpirals = 100;
    int minNGoodSpirals = 2;
    int currentMinNPoints = 99; // between minNPoints and maxNPoints based on distance to first point found
    int nPointsUsed; // nPoints actually used in interpolation
    int imin;
    int imax;
    int jmin;
    int jmax; // current box for this spiral
    /** spiral can do bottom of square */
    boolean doBottom;
    /** spiral can do top of square */
    boolean doTop;
    /** spiral can do left of square */
    boolean doLeft;
    /** spiral can do right of square */
    boolean doRight;
    // reinitialized for each interpolation
    int nSpirals = 0;
    int nPoints = 0;
    int nGoodSpirals = 0;
    boolean spiralHasGoodPoints; // this spiral found a good point
    // running sum of weights; used to determine if point weight is less than minimum
    float weightSum = 0.0f;
    float weightMinFactor = 0.001f;
    boolean weightMinimum;
    /** if point weight is less than 1/distance-squared where distance is the maximum distance considered to be interpolating,
     *  then the point is considered extrapolated; otherwise it is interpolated.
     */
    boolean pointInterpolated = true;

    public ArrayList sortedPoints = new ArrayList(maxNPoints); // points used in interpolation sorted by distance
    static final float nullValue = StsParameters.nullValue;
    static final float largeFloat = StsParameters.largeFloat;

    static final float minWtForInterpolated = 1.0f/16.0f;

    /** implementation provides a loop over rows and columns with appropriate logic for interpolating at each point */
    abstract public float interpolatePoint(int i, int j);
    /** implemention adds a point to the interpolation sortedPoints */
    abstract void addPoint(int i, int j);
    /** implementation does actual interpolation from given points and their weights */
    abstract float getZ(ArrayList sortedPoints);
    /** get the max number of points available for interpolation */
    abstract boolean isPointGrid(int row, int col);
    /** returns name of this surface or grid (used in message). */
    abstract String getGridName();
    /** call if you need to set parameters to values other than the default values */
    public void initializeParameters(int minNPoints, int maxNPoints, int minNGoodSpirals, int maxNSpirals)
    {
        this.minNPoints = minNPoints;
        this.maxNPoints = maxNPoints;
        this.minNGoodSpirals = minNGoodSpirals;
        this.maxNSpirals = maxNSpirals;
    }

    /** set the max number of points available for interpolation. Used in isConverged() criteria.  */
    void setNAvailablePoints()
    {
        nAvailablePoints = 0;
        for(int row = rowMin; row <= rowMax; row++)
            for(int col = colMin; col <= colMax; col++)
                if(isPointGrid(row, col))
                {
                    nAvailablePoints++;
                    if(nAvailablePoints > maxNPoints) return;
                }

    }

    protected void initialize(int iCenter, int jCenter, boolean useGradient, boolean debugInterpolate)
    {
        clearInterpolatingPoints();
        this.iCenter = iCenter;
        this.jCenter = jCenter;
        this.useGradient = useGradient;
        
        if(useGradient)
            minNPoints = 10;
        else
            minNPoints = 3;
        this.debugInterpolate = debugInterpolate;
    }

    public void clearInterpolatingPoints()
    {
        sortedPoints.clear();
        nPoints = 0;
    }

    /** currentMinNPoints is the used in the convergence criteria; it is equal to the distance-squared to the first point found
     *  and is called and set when that first point is found.  It cannot exceed the current number of available points which
     *  are the original ungapped points on the grid.
     * @param distSq the distance-squared to the first point found
     */
    protected void setMinNPoints(int distSq)
    {
        currentMinNPoints = StsMath.minMax(distSq, minNPoints, maxNPoints);
        currentMinNPoints = Math.min(currentMinNPoints, nAvailablePoints);
    }

    protected boolean isConverged()
    {
        if(nPoints >= currentMinNPoints)
        {
            isConverged = true;
            return true;
        }

        if(spiralHasGoodPoints) nGoodSpirals++;
        if(nGoodSpirals < minNGoodSpirals) return false;

        if(weightMinimum)
        {
            isConverged = true;
            return true;
        }
        return false;
    }


    public void interpolateSurface()
    {
        for(int i = rowMin; i <= rowMax; i++)
            for(int j = colMin; j <= colMax; j++)
                interpolatePoint(i, j);
    }

    public float interpolatePoint()
    {
        int i, j;

        try
        {
            imin = iCenter;
            imax = iCenter;
            jmin = jCenter;
            jmax = jCenter;

            clearInterpolatingPoints();
            nSpirals = 0;
            nGoodSpirals = 0;
			weightSum = 0.0f;
            weightMinimum = false;
            isConverged = false;

            if(debugInterpolate)
                StsException.systemDebug(this, "interpolatePoint", "row " + iCenter + " col " + jCenter);  
            while(true)
            {
                spiralHasGoodPoints = false;
                boolean terminate = true;

                doBottom = imin > rowMin;
                if(doBottom) { imin--; terminate = false; }
                doTop = imax < rowMax;
                if(doTop) { imax++; terminate = false; }
                doLeft = jmin > colMin;
                if(doLeft) { jmin--;  terminate = false; }
                doRight = jmax < colMax;
                if(doRight) { jmax++;  terminate = false; }

                if(terminate) break;

                computeNextSpiralWeights(iCenter, jCenter, nSpirals);

                if(doBottom)
                {
                    int max;
                    if(doRight)
                        max = jmax-1;
                    else
                        max = colMax;
                    for(j = jmin; j <= max; j++)
                        addPoint(imin, j);
                }

                if(doRight)
                {
                    int max;
                    if(doTop)
                        max = imax-1;
                    else
                        max = rowMax;
                    for(i = imin; i <= max; i++)
                        addPoint(i, jmax);
                }
                if(doTop)
                {
                    int min;
                    if(doLeft)
                        min = jmin+1;
                    else
                        min = colMin;
                    for(j = min; j <= jmax; j++)
                        addPoint(imax, j);
                }
                if(doLeft)
                {
                    int min;
                    if(doBottom)
                        min = imin+1;
                    else
                        min = rowMin;
                    for(i = min; i <= imax; i++)
                        addPoint(i, jmin);
                }
                if(isConverged()) break;
            }

            nPoints = sortedPoints.size();
            if(nPoints < minNPoints) errorOutput();
            if(nPoints == 0) return nullValue;
            Collections.sort(sortedPoints); // sort points by decreasing weight
            Point nearestPoint = (Point)sortedPoints.get(0);
            pointInterpolated = nearestPoint.wt >= minWtForInterpolated;
            return getZ(sortedPoints);
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "interpolatePoint",
                " for row: " + iCenter + " of " + rowMin + "-" + rowMax +
                " col: " + jCenter + " of " + colMin + "-" + colMax, e);
            return nullValue;
        }
    }

    /** override in subclass if weights are adjusted for depending on ij location and spiral number */
    protected void computeNextSpiralWeights(int i, int j, int nSpiral)
    {
    }

    public boolean getIsConverged() { return isConverged; }

    public void errorOutput()
    {
        StsMessageFiles.errorMessage("Couldn't interpolate point " +
                " for row " + iCenter + " col " + jCenter + " needed: " + minNPoints + " points, found " + nPoints +
                " . See log file.");

        nPoints = sortedPoints.size();
        for(int n = 0; n < nPoints; n++)
        {
            Point point = (Point)sortedPoints.get(n);
            int row = point.row;
            int col = point.col;
            float z = point.z;

            StsMessageFiles.errorMessage(" point " + n + " row: " + row + " col: " + col + " z: " + z);
        }
    }

    static public ArrayList getPoints()
    {
        if(interpolator == null) return null;
        else return interpolator.sortedPoints;
    }

    public class Point implements Comparable
    {
        public int row;
        public int col;
        public float z;
        public float dx;
        public float dy;
        public float wt;

        public Point()
        {
        }

        public Point(int row, int col, float z, float dx, float dy, float wt)
        {
            this.row = row;
            this.col = col;
            this.z = z;
            this.dx = dx;
            this.dy = dy;
            this.wt = wt;
        }

        public int compareTo(Object other)
        {
            float otherWt = ((StsBlockGridInterpolation.Point)other).wt;
            if(wt > otherWt) return -1;
            else if(wt < otherWt) return 1;
            else return 0;
        }

		public String toString()
		{
			return new String("row: " + row + " col: " + col + " z: " + z + " dx: " + dx + " dy: " + dy + " wt: " + wt);
        }
    }
}
