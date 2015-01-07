
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.Utilities.*;

import java.awt.geom.*;

/** a closed polygon used to describe a fault heave on a surface.  It is
    assumed that the points describe a polygon without any internal lines
    intersecting others within in it.  X-Y values are absolute coordinates.
*/
public class StsFaultPolygon extends StsSerialize implements StsSerializable
{
    private double[] x;
    private double[] y;

    private int[] rows;
    private int[] cols;
    private byte[] nullTypes;

    transient private int nPoints;
    transient private Rectangle2D.Double boundingBox;
    transient private Line2D.Double[] segments;
    transient private double yDelta;  // small y shift

    public StsFaultPolygon() { }

	public boolean initialize()
    {
        setNPoints();
        buildBoundingBox();
        buildSegments();
        return true;
    }

    public StsFaultPolygon(double[] x, double[] y)
    {
        // error checks
        if (x == null || y == null || x.length != y.length) return;
        this.x = x;
        this.y = y;
        setNPoints();

        // build a bounding box
        buildBoundingBox();

        // build and save the polygon's line segments
        buildSegments();
    }

    /** see if a point is inside the polygon (uses absolute coordinates) */
    public boolean contains(double xp, double yp)
    {
        if (segments == null || segments.length < 3) return false;
        Point2D.Double p = new Point2D.Double((double)xp, (double)yp);

        // create a reference segment: horizontal line (xp,yp)-(infinity,yp)
        double xpShifted = boundingBox.x - boundingBox.width;
        final Line2D.Double refSegment = new Line2D.Double(xp, yp, xpShifted, yp);

        // if we don't intersect the bounding box, we're done
        if (!refSegment.intersects(boundingBox)) return false;

        // create another reference segment slightly off from horizontal
        final Line2D.Double refSegment2 = new Line2D.Double(xp, yp, xpShifted,
                yp + yDelta);

        // count intersections of the reference segments with the polygon
        int nIntersections = 0;
        for (int i=0; i<nPoints; i++)
        {
            if (!refSegment.intersectsLine(segments[i])) continue;
            if (onSegment(segments[i], p)) return true;
            if (refSegment2.intersectsLine(segments[i])) nIntersections++;
        }

        // odd number of crossings  : the point is inside the polygon
        // even number of crossings : the point is outside the polygon
        return (nIntersections%2 == 1);
    }

    // see if a line contains a point (with floating error StsParameters.roundOff)
    private boolean onSegment(Line2D.Double segment, Point2D point)
    {
        double dist = segment.ptSegDist(point);
        return (dist < StsParameters.roundOff);
    }

    /** accessors */
    public void setX(double[] x) { this.x = x; }
    public double[] getX() { return x; }
    public void setY(double[] y) { this.y = y; }
    public double[] getY() { return y; }
    public void setNPoints(int nPoints) { this.nPoints = nPoints; }
    public int getNPoints() { return nPoints; }
    public void setRows(int[] rows) { this.rows = rows; }
    public int[] getRows() { return rows; }
    public void setCols(int [] cols) { this.cols = cols; }
    public int[] getCols() { return cols; }
    public void setNullTypes(byte[] nullTypes) { this.nullTypes = nullTypes; }
    public byte[] getNullTypes() { return nullTypes; }

    // set no. of points
    private void setNPoints()
    {
        // discard repeated first/last point (if found)
        nPoints = x.length;
        if (x[0] == x[nPoints-1] && y[0] == y[nPoints-1]) nPoints--;
    }

    // build the bounding box
    private void buildBoundingBox()
    {
        if (x == null || y == null) return;
        double xMin = x[0];
        double xMax = x[0];
        double yMin = y[0];
        double yMax = y[0];
        for (int i=1; i<nPoints; i++)
        {
            if (x[i] < xMin) xMin = x[i];
            else if (x[i] > xMax) xMax = x[i];
            if (y[i] < yMin) yMin = y[i];
            else if (y[i] > yMax) yMax = y[i];
        }
        double xWidth = xMax - xMin;
        double yHeight = yMax - yMin;
        boundingBox = new Rectangle2D.Double(xMin, yMin, xWidth, yHeight);
        yDelta = boundingBox.width * StsParameters.roundOff;
    }

    // build the line segments
    private void buildSegments()
    {
        segments = new Line2D.Double[nPoints];
        for (int i=0; i<nPoints; i++)
        {
            int nextI = (i+1) % nPoints;
            segments[i] = new Line2D.Double(x[i], y[i], x[nextI], y[nextI]);
        }
    }


    /** test program */
    public static void main(String[] args)
    {
        final double[] x = { 100.0, 100.0, 200.0, 250.0, 300.0, 300.0, 100.0 };
        final double[] y = { 100.0, 200.0, 200.0, 100.0, 300.0, 0.0,   100.0 };

        StsFaultPolygon fp = new StsFaultPolygon(x, y);
        boolean inside = fp.contains(200.0, 100.0);      // true - inside
        boolean inside2 = fp.contains(200.0, 200.0);     // true - exact point
        boolean inside3 = fp.contains(300.0, 150.0);     // true - on a side
        boolean inside3a = fp.contains(297.0, 150.0);    // true - just inside a side
        boolean inside3b = fp.contains(303.0, 150.0);    // false - just outside a side
        boolean inside4 = fp.contains(400.0, 200.0);     // false
        boolean inside5 = fp.contains(500.0, 500.0);     // false
        boolean inside6 = fp.contains(300.0, 400.0);     // false
        boolean inside7 = fp.contains(400.0, 100.0);     // false
    }
}
