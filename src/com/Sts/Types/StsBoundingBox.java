//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.io.*;

public class StsBoundingBox extends StsMainTimeObject implements Cloneable, Serializable
{
    public float xMin = largeFloat;
    public float yMin = largeFloat;
    public float zMin = largeFloat;
    public float xMax = -largeFloat;
    public float yMax = -largeFloat;
    public float zMax = -largeFloat;

    public double xOrigin, yOrigin;
    public boolean originSet = false;

    transient public boolean initializedXY = false;
    transient public boolean initializedZ = false;

    transient static public final float nullValue = StsParameters.nullValue;
    transient static public final float largeFloat = StsParameters.largeFloat;

    static final long serialVersionUID = -5962170170363197056L;

    static public final int XMIN = 0;
    static public final int XMAX = 1;
    static public final int YMIN = 2;
    static public final int YMAX = 3;
    static public final int ZMIN = 4;
    static public final int ZMAX = 5;
    static public final int CENTER = -1;
    static public final int NONE = -99;

    public StsBoundingBox()
    {
    }

    /** Rotated Bounding Box constructTraceAnalyzer allowing non-persistent construction */
	public StsBoundingBox(boolean persistent)
	{
		super(persistent);
    }

	public StsBoundingBox(boolean persistent, String name)
    {
        super(persistent, name);
    }

    public StsBoundingBox(float xMin, float xMax, float yMin, float yMax)
    {
        xOrigin = 0.0;
        yOrigin = 0.0;
        this.xMin = xMin;
        this.xMax = xMax;
        this.yMin = yMin;
        this.yMax = yMax;
        initializedXY = true;
    }

    public StsBoundingBox(StsXYSurfaceGridable grid)
    {
        xOrigin = grid.getXOrigin();
        yOrigin = grid.getYOrigin();
        xMin = 0.0f;
        yMin = 0.0f;
        xMax = grid.getXSize(); // only used if unrotated
        yMax = grid.getYSize(); // only used if unrotated
        zMin = grid.getZMin();
        zMax = grid.getZMax();
        initializedXY = true;
        initializedZ = true;
    }

    public StsBoundingBox(double xOrigin, double yOrigin)
    {
        this.xOrigin = xOrigin;
        this.yOrigin = yOrigin;
    }

    public StsBoundingBox(StsPoint[] points, double xOrigin, double yOrigin)
    {
        addPoints(points, xOrigin, yOrigin);
    }


    public StsBoundingBox(double[][] points)
    {
		for(double[] point : points)
        	addPoint(point);
    }

    public StsBoundingBox(StsPoint[] points)
    {
        addPoints(points);
    }

    public StsBoundingBox getClone()
    {
        try
        {
            return (StsBoundingBox)this.clone();
        }
        catch (Exception e)
        {
            StsException.outputException("StsBoundingBox.getClone(grid) failed.", e, StsException.WARNING);
            return null;
        }
    }

    public void initialize(float xMin, float xMax, float yMin, float yMax, float zMin, float zMax)
    {
        initialize(xMin, xMax, yMin, yMax, zMin, zMax, true);
    }

    public void initialize(float xMin, float xMax, float yMin, float yMax, float zMin, float zMax, boolean initialized)
    {
        this.xMin = xMin;
        this.xMax = xMax;
        this.yMin = yMin;
        this.yMax = yMax;
        this.zMin = zMin;
        this.zMax = zMax;
        this.initializedXY = initialized;
        this.initializedZ = initialized;
    }

    public boolean initialize(StsModel model)
    {
        initializedXY = true;
        initializedZ = true;
        return true;
    }

    /** classInitialize the boundingBox with a new origin and reset the limits.
     *  This occurs when the first object generated boundingBox is added to the project.
     */
    public void initializeOriginAndRange(double xOrigin, double yOrigin)
    {
        this.xOrigin = xOrigin;
        this.yOrigin = yOrigin;
        originSet = true;
        resetRange();
    }

	public void reinitializeBoundingBox()
	{
		resetRange();
	}

    public void resetRange()
    {
       resetXYRange();
       resetZRange();
    }

    public void resetXYRange()
    {
        xMin = largeFloat;
        yMin = largeFloat;
        xMax = -largeFloat;
        yMax = -largeFloat;
        initializedXY = false;
    }

    public void resetZRange()
    {
        zMin = largeFloat;
        zMax = -largeFloat;
        initializedZ = false;
    }

    public boolean isInsideXY(float[] xy)
    {
        return isInsideXY(xy[0], xy[1]);
    }

    public boolean isInsideXY(float x, float y)
    {
        return x >= xMin && x <= xMax && y >= yMin && y <= yMax;
    }

    public boolean isInsideXY(StsPoint p)
    {
        return isInsideXY(p.v[0], p.v[1]);
    }

    public boolean isInsideXYZ(StsPoint p)
    {
        return isInsideXYZ(p.v[0], p.v[1], p.v[2]);
    }

    public boolean isInsideXYZ(float x, float y, float z)
    {
        if (x < xMin) return false;
        if (x > xMax) return false;
        if (y < yMin) return false;
        if (y > yMax) return false;
        if (z < zMin) return false;
        if (z > zMax) return false;
        return true;
   }
   public boolean isInsideXYZ(float[] xyz)
   {
	   return isInsideXYZ(xyz[0], xyz[1], xyz[2]);
    }

    public float getXSize()
    {
		return xMax - xMin;
    }

    public float getYSize()
    {
		return yMax - yMin;
    }

    public float[] getRelativeXY(double x, double y)
    {
        return new float[]
            { (float) (x - xOrigin), (float) (y - yOrigin)};
    }
/*
    public void adjustRange(float[] range)
    {
        xMin = Math.min(xMin, range[0]);
        xMax = Math.max(xMax, range[1]);
        yMin = Math.min(yMin, range[2]);
        yMax = Math.max(yMax, range[3]);
        zMin = Math.min(zMin, range[4]);
        zMax = Math.max(zMax, range[5]);
        initialized = true;
        initializedZ = true;
    }
*/
   public void adjustXYZPosition(StsPoint dPoint)
    {
        float[] dxyz = dPoint.getXYZorT();
        xMin += dxyz[0];
        xMax += dxyz[0];
        yMin += dxyz[1];
        yMax += dxyz[1];
        zMin += dxyz[2];
        zMax += dxyz[2];
    }

	public void addPoint(double[] xyz)
	{
		addPoint(new float[] { (float)xyz[0], (float)xyz[1], (float)xyz[2] } );
	}

    public void addPoint(float[] xyz)
    {
        if(xyz == null) return;

        if (xyz.length > 2)
        {
            float z = xyz[2];
            if (z == nullValue)
            {
                return;
            }
            zMin = Math.min(z, zMin);
            zMax = Math.max(z, zMax);
            initializedZ = true;
        }
        float x = xyz[0];
        float y = xyz[1];

        xMin = Math.min(x, xMin);
        xMax = Math.max(x, xMax);
        yMin = Math.min(y, yMin);
        yMax = Math.max(y, yMax);
        initializedXY = true;
    }

    public void addPoint(float[] xyz, float dXOrigin, float dYOrigin)
    {
		float x = xyz[0];
		float y = xyz[1];
		float z = xyz[2];
		addPoint(x, y, z, dXOrigin, dYOrigin);
    }

	public void addPoint(float x, float y, float z, float dXOrigin, float dYOrigin)
	{
		if (z == nullValue) return;

        x += dXOrigin;
        y += dYOrigin;

        xMin = Math.min(x, xMin);
        xMax = Math.max(x, xMax);
        yMin = Math.min(y, yMin);
        yMax = Math.max(y, yMax);
        initializedXY = true;
        zMin = Math.min(z, zMin);
        zMax = Math.max(z, zMax);
        initializedZ = true;
    }

    public void addPoint(StsGridSectionPoint point)
    {
        addPoint(point.getXYZorT());
    }

    public void addPoints(StsPoint[] points)
    {
        if (points == null)
        {
            return;
        }
        for (int n = 0; n < points.length; n++)
        {
            StsPoint point = points[n];
            if (point == null)
            {
                return;
            }
            addPoint(point.v);
        }
    }

    public void addPoints(double[][] points)
    {
        int nPoints = points.length;
        if (nPoints == 0)
        {
            return;
        }
        int nCoors = Math.min(points[0].length, 3);
        if (nCoors == 0)
        {
            return;
        }

        for (int n = 0; n < nCoors; n++)
        {
            float min = StsParameters.largeFloat;
            float max = -StsParameters.largeFloat;

            for (int p = 0; p < nPoints; p++)
            {
                float value = (float) points[p][n];
                min = Math.min(min, value);
                max = Math.max(max, value);
            }

            switch (n)
            {
                case 0:
                    xMin = min;
                    xMax = max;
                    break;
                case 1:
                    yMin = min;
                    yMax = max;
                    break;
                case 2:
                    zMin = min;
                    zMax = max;
            }
        }
        initializedXY = true;
        if(nCoors > 2) initializedZ = true;
    }

    /** This point is offset from its own origin.
     *  Add it to this boundingBox which has its own origin
     */

    public void addPoint(StsPoint point, double xOrigin, double yOrigin)
    {
        if (point == null)
        {
            return;
        }
        float dXOrigin = (float) (xOrigin - this.xOrigin);
        float dYOrigin = (float) (yOrigin - this.yOrigin);
        addPoint(point.v, dXOrigin, dYOrigin);
    }

    public void addXY(StsPoint point, double xOrigin, double yOrigin)
    {
        if (point == null) return;
        float dXOrigin = (float) (xOrigin - this.xOrigin);
        float dYOrigin = (float) (yOrigin - this.yOrigin);
        float x = point.v[0] + dXOrigin;
        float y = point.v[1] + dYOrigin;

        xMin = Math.min(x, xMin);
        xMax = Math.max(x, xMax);
        yMin = Math.min(y, yMin);
        yMax = Math.max(y, yMax);
        initializedXY = true;
    }


    public void addXY(double x, double y, double xOrigin, double yOrigin)
    {
        float dXOrigin = (float) (xOrigin - this.xOrigin);
        float dYOrigin = (float) (yOrigin - this.yOrigin);
        x += dXOrigin;
        y += dYOrigin;

        xMin = (float)Math.min(x, xMin);
        xMax = (float)Math.max(x, xMax);
        yMin = (float)Math.min(y, yMin);
        yMax = (float)Math.max(y, yMax);
        initializedXY = true;
    }

    private boolean sanityCheck(double x, double y, double limit)
    {
        if(xMin == StsParameters.largeFloat)     // First point in project always passes sanity check
            return true;
        // check point against current bounding box
        double minDx = Math.abs(x - xMin);
        double maxDx = Math.abs(xMax - x);
        double minDy = Math.abs(y - yMin);
        double maxDy = Math.abs(yMax - y);
        if((minDx > limit) || (maxDx > limit) || (minDy > limit) || (maxDy > limit))
            return false;
        return true;
    }
    public boolean sanityCheck(double x, double y, double xOrigin, double yOrigin, double limit)
    {
        float dXOrigin = (float) (xOrigin - this.xOrigin);
        float dYOrigin = (float) (yOrigin - this.yOrigin);
        x += dXOrigin;
        y += dYOrigin;

        // check point against current bounding box
        return sanityCheck(x,y, limit);
    }

    public boolean sanityCheck(StsPoint point, double xOrigin, double yOrigin, double limit)
    {
        if (point == null)
            return false;
        float dXOrigin = (float) (xOrigin - this.xOrigin);
        float dYOrigin = (float) (yOrigin - this.yOrigin);
        float x = point.v[0] + dXOrigin;
        float y = point.v[1] + dYOrigin;

        // check point against current bounding box
        return sanityCheck(x,y, limit);
    }

    public boolean sanityCheck(StsBoundingBox box, double limit)
    {
        float dx = (float) (box.xOrigin - xOrigin);
        float dy = (float) (box.yOrigin - yOrigin);
        double x = dx + box.xMin;
        double y = dy + box.yMin;

        // check point against current bounding box
        return sanityCheck(x,y, limit);
    }

    public void addPoint(float x, float y, float z, double xOrigin, double yOrigin)
	 {
		 float dXOrigin = (float) (xOrigin - this.xOrigin);
		 float dYOrigin = (float) (yOrigin - this.yOrigin);
		 addPoint(x, y, z, dXOrigin, dYOrigin);
	 }

    public void addPoints(StsPoint[] points, double xOrigin, double yOrigin)
    {
        if (points == null)
        {
            return;
        }

        float dXOrigin = (float) (xOrigin - this.xOrigin);
        float dYOrigin = (float) (yOrigin - this.yOrigin);
        for (int n = 0; n < points.length; n++)
        {
            StsPoint point = points[n];
            if (point == null)
            {
                return;
            }
            addPoint(point.v, dXOrigin, dYOrigin);
        }
    }

    /** reset this boundingBox to include box.
     *  This method should not be used for a rotated bounding box!
     */
    public void addUnrotatedBoundingBox(StsBoundingBox box)
    {
        if (!originSet)
        {
            initializeOriginAndRange(box.xOrigin, box.yOrigin);
            if (initializedXY)
            {
                xMin = Math.min(xMin, box.xMin);
                xMax = Math.max(xMax, box.xMax);
                yMin = Math.min(yMin, box.yMin);
                yMax = Math.max(yMax, box.yMax);
            }
            else
            {
                xMin = box.xMin;
                xMax = box.xMax;
                yMin = box.yMin;
                yMax = box.yMax;
                initializedXY = true;
            }
        }
        else
        {
            float dx = (float) (box.xOrigin - xOrigin);
            float dy = (float) (box.yOrigin - yOrigin);
            box.xOrigin = xOrigin;
            box.yOrigin = yOrigin;
            box.xMin += dx;
            box.xMax += dx;
            box.yMin += dy;
            box.yMax += dy;
            if (initializedXY)
            {
                xMin = Math.min(xMin, box.xMin);
                xMax = Math.max(xMax, box.xMax);
                yMin = Math.min(yMin, box.yMin);
                yMax = Math.max(yMax, box.yMax);
            }
            else
            {
                xMin = box.xMin;
                xMax = box.xMax;
                yMin = box.yMin;
                yMax = box.yMax;
                initializedXY = true;
            }
        }
        if (initializedZ)
        {
            zMin = Math.min(zMin, box.zMin);
            zMax = Math.max(zMax, box.zMax);
        }
        else
        {
            zMin = box.zMin;
            zMax = box.zMax;
            initializedZ = true;
        }
    }

    public void intersectBoundingBox(StsBoundingBox box)
    {
        if (!originSet)
        {
            // Don't want to reset the ranges as well or else zmin and zmax will always be 1E-10 and -1E-10
            // because of these statements below; zMin = Math.max(zMin, box.zMin); zMax = Math.min(zMax, box.zMax);
            // initializeOriginAndRange(box.xOrigin, box.yOrigin);
            originSet = true;
            if(initializedXY)
            {
                xMin = Math.max(xMin, box.xMin);
                xMax = Math.min(xMax, box.xMax);
                yMin = Math.max(yMin, box.yMin);
                yMax = Math.min(yMax, box.yMax);
            }
            else
            {
                xMin = box.xMin;
                xMax = box.xMax;
                yMin = box.yMin;
                yMax = box.yMax;
                initializedXY = true;
            }
        }
        else
        {
            float dx = (float) (box.xOrigin - xOrigin);
            float dy = (float) (box.yOrigin - yOrigin);
            box.xOrigin = xOrigin;
            box.yOrigin = yOrigin;
            box.xMin += dx;
            box.xMax += dx;
            box.yMin += dy;
            box.yMax += dy;
            if (initializedXY)
            {
                xMin = Math.max(xMin, box.xMin);
                xMax = Math.min(xMax, box.xMax);
                yMin = Math.max(yMin, box.yMin);
                yMax = Math.min(yMax, box.yMax);
            }
            else
            {
                xMin = box.xMin;
                xMax = box.xMax;
                yMin = box.yMin;
                yMax = box.yMax;
                initializedXY = true;
            }
        }
        if(initializedZ)
        {
            zMin = Math.max(zMin, box.zMin);
            zMax = Math.min(zMax, box.zMax);
        }
        else
        {
            zMin = box.zMin;
            zMax = box.zMax;
            initializedZ = true;
        }
    }
    public boolean intersectsBoundingBox(StsBoundingBox box)
    {
		double xMinI = Math.max(box.getAbsXMin(), getAbsXMin());
		double xMaxI = Math.min(box.getAbsXMax(), getAbsXMax());
		if(xMaxI <= xMinI) return false;
		double yMinI = Math.max(box.getAbsYMin(), getAbsYMin());
		double yMaxI = Math.min(box.getAbsYMax(), getAbsYMax());
		if(yMaxI <= yMinI) return false;
		float zMinI = Math.max(box.zMin, zMin);
		float zMaxI = Math.min(box.zMax, zMax);
		return zMaxI > zMinI;
    }

    public void initialize(StsBoundingBox box)
    {
//		this.xOrigin = box.xOrigin;
//		this.yOrigin = box.yOrigin;
//		this.originSet = true;
        xMin = box.xMin;
        xMax = box.xMax;
        yMin = box.yMin;
        yMax = box.yMax;
        zMin = box.zMin;
        zMax = box.zMax;
        initializedXY = true;
        initializedZ = true;
    }

    public boolean clipLine(StsPoint p0, StsPoint p1)
    {
		double[] xyz0 = new double[] { p0.v[0], p0.v[1], p0.v[2] };
		double[] xyz1 = new double[] { p1.v[0], p1.v[1], p1.v[2] };
		return clipLine(xyz0, xyz1);
    }

	public boolean clipLine(double[] xyz0, double[] xyz1)
	{
        double max, min;

        /* clip ends between xMin and xMax 				*/

        if (xyz0[0] < xyz1[0])
        {
            clipXLine(xyz0, xyz1);
        }
        else
        {
            clipXLine(xyz1, xyz0);

            /* check if clipped line is below yMin or above yMax	*/

        }
        max = Math.max(xyz0[1], xyz1[1]);
        min = Math.min(xyz0[1], xyz1[1]);

        if (max <= yMin || min >= yMax)
        {
            return false;
        }

        /* clip ends between yMin and yMax 					*/

        if (xyz0[1] < xyz1[1])
        {
            clipYLine(xyz0, xyz1);
        }
        else
        {
            clipYLine(xyz1, xyz0);

            /* check if clipped line is below xMin or above xMax	*/

        }
        max = Math.max(xyz0[0], xyz1[0]);
        min = Math.min(xyz0[0], xyz1[0]);

        if (max <= xMin || min >= xMax)
        {
            return false;
        }
        else
        {
            return true;
        }
    }

    // Line is increasing in X-index from start to end and extends beyond
    // an X-direction boundary: clip it.

    private void clipXLine(StsPoint start, StsPoint end)
    {
        float f;

        if (start.v[0] < xMin && end.v[0] > xMin)
        {
            f = (xMin - start.v[0]) / (end.v[0] - start.v[0]);
            start.v[1] = start.v[1] + f * (end.v[1] - start.v[1]);
            start.v[0] = xMin;
        }

        if (end.v[0] > xMax && start.v[0] < xMax)
        {
            f = (xMax - start.v[0]) / (end.v[0] - start.v[0]);
            end.v[1] = start.v[1] + f * (end.v[1] - start.v[1]);
            end.v[0] = xMax;
        }
    }

    private void clipXLine(double[] start, double[] end)
    {
        double f;

        if (start[0] < xMin && end[0] > xMin)
        {
            f = (xMin - start[0]) / (end[0] - start[0]);
            start[0] = xMin;
            start[1] = start[1] + f * (end[1] - start[1]);
            start[2] = start[2] + f * (end[2] - start[2]);
        }

        if (end[0] > xMax && start[0] < xMax)
        {
            f = (xMax - start[0]) / (end[0] - start[0]);
            end[0] = xMax;
            end[1] = start[1] + f * (end[1] - start[1]);
            end[2] = start[2] + f * (end[2] - start[2]);
        }
    }

    // Line is increasing in Y-index from start to end and extends beyond
    // a Y-direction boundary: clip it.

    private void clipYLine(StsPoint start, StsPoint end)
    {
        float f;

        if (start.v[1] < yMin && end.v[1] > yMin)
        {
            f = (yMin - start.v[1]) / (end.v[1] - start.v[1]);
            start.v[0] = start.v[0] + f * (end.v[0] - start.v[0]);
            start.v[1] = yMin;
        }

        if (end.v[1] > yMax && start.v[1] < yMax)
        {
            f = (yMax - start.v[1]) / (end.v[1] - start.v[1]);
            end.v[0] = start.v[0] + f * (end.v[0] - start.v[0]);
            end.v[1] = yMax;
        }
    }

    private void clipYLine(double[] start, double[] end)
    {
        double f;

        if (start[1] < yMin && end[1] > yMin)
        {
            f = (yMin - start[1]) / (end[1] - start[1]);
            start[1] = yMin;
            start[0] = start[0] + f * (end[0] - start[0]);
            start[2] = start[2] + f * (end[2] - start[2]);
        }

        if (end[1] > yMax && start[1] < yMax)
        {
            f = (yMax - start[1]) / (end[1] - start[1]);
            end[1] = yMax;
            end[0] = start[0] + f * (end[0] - start[0]);
            end[2] = start[2] + f * (end[2] - start[2]);
        }
    }

    public void addRotatedBoundingBox(StsRotatedBoundingBox boundingBox)
    {
        double[][] corners = new double[4][];
        corners[0] = boundingBox.getAbsoluteXY(boundingBox.xMin, boundingBox.yMin);
        corners[1] = boundingBox.getAbsoluteXY(boundingBox.xMax, boundingBox.yMin);
        corners[2] = boundingBox.getAbsoluteXY(boundingBox.xMax, boundingBox.yMax);
        corners[3] = boundingBox.getAbsoluteXY(boundingBox.xMin, boundingBox.yMax);

        for (int n = 0; n < 4; n++)
        {
            float x = (float) (corners[n][0] - xOrigin);
            float y = (float) (corners[n][1] - yOrigin);
            xMin = Math.min(xMin, x);
            xMax = Math.max(xMax, x);
            yMin = Math.min(yMin, y);
            yMax = Math.max(yMax, y);
        }
        zMin = Math.min(zMin, boundingBox.zMin);
        zMax = Math.max(zMax, boundingBox.zMax);
    }

    public String getLabel()
    {
        return toDetailString();
    }

    public String toDetailString()
    {
        return super.toString() + " xMin: " + xMin + " xMax: " + xMax + " yMin: " + yMin + " yMax: " + yMax + " zMin: " + zMin + " zMax: " + zMax;
    }

    public static void main(String[] args)
    {
        float xMin = -1.0f;
        float xMax = 1.0f;
        float yMin = -1.0f;
        float yMax = 1.0f;

        StsBoundingBox box = new StsBoundingBox( -1.0f, 1.0f, -1.0f, 1.0f);
        System.out.println("xMin: " + xMin + " xMax: " + xMax +
                           " yMin: " + yMin + " yMax: " + yMax);

        StsPoint p0 = new StsPoint( -2.0f, -2.0f, -2.0f);
        StsPoint p1 = new StsPoint(2.0f, 2.0f, 2.0f);
        System.out.println("input points: " + p0.toString() + " " + p1.toString());

        box.clipLine(p0, p1);
        System.out.println("output points: " + p0.toString() + " " + p1.toString());
    }

    public void setXOrigin(double xOrigin) { this.xOrigin = xOrigin; }
    public void setYOrigin(double yOrigin) { this.yOrigin = yOrigin; }

    public boolean setOrigin(double xOrigin, double yOrigin)
    {
        if (originSet) return false;
        initializeOriginAndRange(xOrigin, yOrigin);
        return true;
    }

    public void setZRange(float zMin, float zMax)
    {
        this.zMin = zMin;
        this.zMax = zMax;
    }

	public boolean sameAs(StsBoundingBox otherBox)
	{
		return sameAs(otherBox, true);
	}

   public boolean sameAs(StsBoundingBox otherBox, boolean checkZ)
    {
        if (!StsMath.sameAs(xOrigin, otherBox.xOrigin))
        {
            return false;
        }
        if (!StsMath.sameAs(yOrigin, otherBox.yOrigin))
        {
            return false;
        }
        if (!StsMath.sameAs(xMin, otherBox.xMin))
        {
            return false;
        }
        if (!StsMath.sameAs(xMax, otherBox.xMax))
        {
            return false;
        }
        if (!StsMath.sameAs(yMin, otherBox.yMin))
        {
            return false;
        }
        if (!StsMath.sameAs(yMax, otherBox.yMax))
        {
            return false;
        }
        if (checkZ && !StsMath.sameAs(zMin, otherBox.zMin))
        {
            return false;
        }
        if (checkZ && !StsMath.sameAs(zMax, otherBox.zMax))
        {
            return false;
        }
        return true;
    }

    /** Get the x center of the box */
    public float getXCenter()
    {return (xMin + xMax) / 2;
    }

    /** Get the y center of the box */
    public float getYCenter()
    {return (yMin + yMax) / 2;
    }

    /** Get the maximum Project dimensions
     * @returns maximum projection distance in X, Y or Z
     */
    public float getDimensions()
    {
        return StsMath.max3( (xMax - xMin), (yMax - yMin), (zMax - zMin));
    }

    public double getXOrigin()
    {return xOrigin;
    }

    public double getYOrigin()
    {return yOrigin;
    }

    public void setZMin(float zMin)
    {this.zMin = zMin;
    }

    public void setZMax(float zMax) { this.zMax = zMax; }

    public float getZMin()
    {
        return zMin;
    }

    public float getZMin(int row)
    {
        return zMin;
    }

    public float getZMax()
    {
        return zMax;
    }

    public float getYMin()
    {
        return yMin;
    }

    public float getYMax()
    {
        return yMax;
    }

    public void setYMax(float yMax)
    {
        this.yMax = yMax;
    }

    public void setYMin(float yMin)
    {
        this.yMin = yMin;
    }

    public void setXMin(float xMin)
    {
		this.xMin = xMin;
    }

    public void setXMax(float xMax)
    {
        this.xMax = xMax;
    }

    public float getXMin()
    {
        return xMin;
    }

    public float getXMax()
    {
        return xMax;
    }

	public double getAbsXMin() { return xOrigin + xMin; }
	public double getAbsXMax() { return xOrigin + xMax; }
	public double getAbsYMin() { return yOrigin + yMin; }
	public double getAbsYMax() { return yOrigin + yMax; }
	public void setXMin(Float xMin) { setXMin(xMin.floatValue()); }
	public void setXMax(Float xMax) { setXMax(xMax.floatValue()); }
	public void setYMin(Float yMin) { setYMin(yMin.floatValue()); }
	public void setYMax(Float yMax) { setYMax(yMax.floatValue()); }
	public void setZMin(Float zMin) { setZMin(zMin.floatValue()); }
	public void setZMax(Float zMax) { setZMax(zMax.floatValue()); }

    public void displayBoundingBox(GL gl, StsColor stsColor, float lineWidth)
    {
        float[] vec = new float[3];

        stsColor.setGLColor(gl);

        gl.glLineWidth(lineWidth);

        // Draw the bottom outer boundaries
        gl.glBegin(GL.GL_LINE_LOOP);
        {
            vec[2] = zMax;
            vec[0] = xMin;
            vec[1] = yMin;
            gl.glVertex3fv(vec, 0);
            vec[0] = xMin;
            vec[1] = yMax;
            gl.glVertex3fv(vec, 0);
            vec[0] = xMax;
            vec[1] = yMax;
            gl.glVertex3fv(vec, 0);
            vec[0] = xMax;
            vec[1] = yMin;
            gl.glVertex3fv(vec, 0);
        }
        gl.glEnd();

        // Draw the top outer boundaries
        gl.glBegin(GL.GL_LINE_LOOP);
        {
            vec[2] = zMin;
            vec[0] = xMin;
            vec[1] = yMin;
            gl.glVertex3fv(vec, 0);
            vec[0] = xMin;
            vec[1] = yMax;
            gl.glVertex3fv(vec, 0);
            vec[0] = xMax;
            vec[1] = yMax;
            gl.glVertex3fv(vec, 0);
            vec[0] = xMax;
            vec[1] = yMin;
            gl.glVertex3fv(vec, 0);
        }
        gl.glEnd();

        // Draw verticals
        gl.glBegin(GL.GL_LINES);
        {
            vec[0] = xMin;
            vec[1] = yMin;
            vec[2] = zMin;
            gl.glVertex3fv(vec, 0);
            vec[2] = zMax;
            gl.glVertex3fv(vec, 0);

            vec[1] = yMax;
            vec[2] = zMin;
            gl.glVertex3fv(vec, 0);
            vec[2] = zMax;
            gl.glVertex3fv(vec, 0);

            vec[0] = xMax;
            vec[2] = zMin;
            gl.glVertex3fv(vec, 0);
            vec[2] = zMax;
            gl.glVertex3fv(vec, 0);

            vec[1] = yMin;
            vec[2] = zMin;
            gl.glVertex3fv(vec, 0);
            vec[2] = zMax;
            gl.glVertex3fv(vec, 0);
        }
        gl.glEnd();
    }

    public void drawCornerPoints(GL gl, StsColor color, boolean isPicking)
    {
        float[] xyz = new float[3];

        xyz[0] = xMin;
        xyz[1] = yMin;
        xyz[2] = zMin;

        if(isPicking) gl.glPushName(0);
        StsGLDraw.drawPoint(gl, xyz, color, 4);

        xyz[0] = xMax;
        if(isPicking) { gl.glPopName(); gl.glPushName(1); }
        StsGLDraw.drawPoint(gl, xyz, color, 4);

        xyz[1] = yMax;
        if(isPicking) { gl.glPopName(); gl.glPushName(2); }
        StsGLDraw.drawPoint(gl, xyz, color, 4);

        xyz[0] = xMin;
        if(isPicking) { gl.glPopName(); gl.glPushName(3); }
        StsGLDraw.drawPoint(gl, xyz, color, 4);

        xyz[2] = zMax;
        xyz[1] = yMin;

        if(isPicking) gl.glPushName(4);
        StsGLDraw.drawPoint(gl, xyz, color, 4);

        xyz[0] = xMax;
        if(isPicking) { gl.glPopName(); gl.glPushName(5); }
        StsGLDraw.drawPoint(gl, xyz, color, 4);

        xyz[1] = yMax;
        if(isPicking) { gl.glPopName(); gl.glPushName(6); }
        StsGLDraw.drawPoint(gl, xyz, color, 4);

        xyz[0] = xMin;
        if(isPicking) { gl.glPopName(); gl.glPushName(7); }
        StsGLDraw.drawPoint(gl, xyz, color, 4);

        if(isPicking) gl.glPopName();
    }

    public void drawFacePoints(GL gl, boolean isPicking)
    {
        float[] xyz = new float[3];

        xyz[2] = (zMin + zMax)/2;

        xyz[0] = xMin;
        xyz[1] = (yMin+yMax)/2;
        if(isPicking) gl.glPushName(0);
        StsGLDraw.drawPoint(gl, xyz, StsColor.WHITE, 4);

        xyz[0] = xMax;
        if(isPicking) { gl.glPopName(); gl.glPushName(1); }
        StsGLDraw.drawPoint(gl, xyz, StsColor.WHITE, 4);

        xyz[0] = (xMin + xMax)/2;
        xyz[1] = yMin;
        if(isPicking) { gl.glPopName(); gl.glPushName(2); }
        StsGLDraw.drawPoint(gl, xyz, StsColor.WHITE, 4);

        xyz[1] = yMax;
        if(isPicking) { gl.glPopName(); gl.glPushName(3); }
        StsGLDraw.drawPoint(gl, xyz, StsColor.WHITE, 4);

        xyz[1] = (yMin + yMax)/2;
        xyz[2] = zMin;
        if(isPicking) { gl.glPopName(); gl.glPushName(4); }
        StsGLDraw.drawPoint(gl, xyz, StsColor.WHITE, 4);

        xyz[2] = zMax;
        if(isPicking) { gl.glPopName(); gl.glPushName(5); }
        StsGLDraw.drawPoint(gl, xyz, StsColor.WHITE, 4);

        if(isPicking) gl.glPopName();
    }

    public float[] getFaceCenter(int faceIndex)
    {
        switch(faceIndex)
        {
            case XMIN:
                return new float[] { xMin, (yMin + yMax)/2, (zMin + zMax)/2 };
            case XMAX:
                return new float[] { xMax, (yMin + yMax)/2, (zMin + zMax)/2 };
            case YMIN:
                return new float[] { (xMin + xMax)/2, yMin, (zMin + zMax)/2 };
            case YMAX:
               return new float[] { (xMin + xMax)/2, yMax, (zMin + zMax)/2 };
           case ZMIN:
                return new float[] { (xMin + xMax)/2, (yMin + yMax)/2, zMin };
            case ZMAX:
                return new float[] { (xMin + xMax)/2, (yMin + yMax)/2, zMax };
            default:
                return null;
        }
    }

    public void adjustRange(int faceIndex, float adjustment)
    {
        switch(faceIndex)
        {
            case XMIN:
                if((xMin + adjustment) < getXMax())
                    setXMin(xMin += adjustment);
                break;
            case XMAX:
                if((xMax + adjustment) > getXMin())
                    setXMax(xMax += adjustment);
                break;
            case YMIN:
                if((yMin + adjustment) < getYMax())
                    setYMin(yMin += adjustment);
                break;
            case YMAX:
               if((yMax + adjustment) > getYMin())
                   setYMax(yMax += adjustment);
               break;
          case ZMIN:
                if((zMin + adjustment) < getZMax())
                    setZMin(zMin += adjustment);
                break;
            case ZMAX:
                if((zMax + adjustment) > getZMin())
                    setZMax(zMax += adjustment);
                break;
        }
    }

	public double[][] normalizePoints(double[][] points)
	{
		double[] scale = new double[3];
		double[] offset = new double[3];
		scale[0] = 1.0/(xMax - xMin);
		offset[0] = -xMin*scale[0];
		scale[1] = 1.0/(yMax - yMin);
		offset[1] = -yMin*scale[1];
		scale[2] = 1.0/(zMax - zMin);
		offset[2] = -zMin*scale[2];
		double[][] normalizedPoints = StsMath.copy(points);
		int nPoints = points.length;
		for(int n = 0; n < nPoints; n++)
			for(int i = 0; i < 3; i++)
			normalizedPoints[n][i] = offset[i] + scale[i]*points[n][i];
		return normalizedPoints;
	}

	public void addBorderXY(float borderSize)
	{
		xMin -= borderSize;
		xMax += borderSize;
		yMin -= borderSize;
		yMax += borderSize;
	}
}
