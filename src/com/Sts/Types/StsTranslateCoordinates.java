
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.Utilities.*;
/** Defines a coordinate system for a line between two points */

public class StsTranslateCoordinates
{
    /**  	@param a is unit normal along line between first and last points
      *		@param h is unit normal orthogonal to a and horizontal thru first point
      *		@param v is unit normal orthogonal to a and h.
      *     @param o is origin (first point).
      *     @param length is distance between origin and last point
      */

    public StsPoint a, h, v, o;
	public StsPoint f;
    public float length;

    public StsTranslateCoordinates(StsPoint firstPoint, StsPoint lastPoint)
    {
        o = new StsPoint(firstPoint);
        a = StsPoint.subPointsStatic(lastPoint, firstPoint);
        length = a.normalizeReturnLength();
   /*
        if(length <= 0.0f)
            throw new StsException(StsException.WARNING, "Length is zero in LineCoordinates constructTraceAnalyzer");
   */
        h = new StsPoint(a.v[1], -a.v[0], 0.0f);
        h.normalize();
        v = new StsPoint(3);
        v.crossProduct(a, h);
    }

    public StsTranslateCoordinates(StsPoint firstPoint, StsPoint lastPoint, StsPoint point)
    {
        o = new StsPoint(firstPoint);
        a = StsPoint.subPointsStatic(lastPoint, firstPoint);
        length = a.normalizeReturnLength();
        h = new StsPoint(a.v[1], -a.v[0], 0.0f);
        h.normalize();

        v = new StsPoint(3);
        v.crossProduct(a, h);

		StsPoint dPoint = StsPoint.subPointsStatic(point, o);

		f = new StsPoint(3);
        f.v[0] = dPoint.dot(a);
        f.v[1] = dPoint.dot(h);
        f.v[2] = dPoint.dot(v);
    }

    static public void computeTranslatedPoints(StsPoint[] oldPoints, StsTranslateCoordinates oldCoor,
                                               StsPoint[] newPoints, StsTranslateCoordinates newCoor)
        throws StsException
    {
        if(oldPoints == null)
        {
            StsException.systemError("StsTranslateCoordinates.computeTranslatedPoint() failed." +
                " oldPoints null in computeTranslatedPoints");
            return;
        }
        if(oldCoor == null)
        {
            StsException.systemError("StsTranslateCoordinates.computeTranslatedPoint() failed." +
                " oldCoor null in computeTranslatedPoints");
            return;
        }
        if(newPoints == null)
        {
            StsException.systemError("StsTranslateCoordinates.computeTranslatedPoint() failed." +
                " newPoints null in computeTranslatedPoints");
            return;
        }
        if(newCoor == null)
        {
            StsException.systemError("StsTranslateCoordinates.computeTranslatedPoint() failed." +
                " newCoor null in computeTranslatedPoints");
            return;
        }

        float scaleFactor = newCoor.length/oldCoor.length;

        /** now compute new internal points by scaling from old points */
        for(int n = 0; n < oldPoints.length; n++)
        {
            StsPoint dPoint = StsPoint.subPointsStatic(oldPoints[n], oldCoor.o);

            float fa = scaleFactor*dPoint.dot(oldCoor.a, 3);
            float fh = scaleFactor*dPoint.dot(oldCoor.h, 3);
            float fv = scaleFactor*dPoint.dot(oldCoor.v, 3);

            float[] v = new float[4];

            for(int i = 0; i < 3; i++)
                v[i] = fa*newCoor.a.v[i] + fh*newCoor.h.v[i] + fv*newCoor.v.v[i] + newCoor.o.v[i];

            if(newPoints[n] == null)
                newPoints[n] = new StsPoint(v);
            else
                newPoints[n].setValues(v);
        }
    }

    static public StsPoint[] computeTranslatedPoints(StsPoint[] oldPoints, StsTranslateCoordinates oldCoor,
                                                     StsPoint firstPoint, StsPoint lastPoint)
        throws StsException
    {
        StsTranslateCoordinates newCoor = new StsTranslateCoordinates(firstPoint, lastPoint);

        StsPoint[] newPoints = new StsPoint[oldPoints.length];

        computeTranslatedPoints(oldPoints, oldCoor, newPoints, newCoor);

        return newPoints;
    }

	// computes oldCoordinates of oldPoint and returns xyz of point translated to newCoorinates
    static public StsPoint computeTranslatedPoint(StsPoint oldPoint, StsTranslateCoordinates oldCoor,
                                                  StsTranslateCoordinates newCoor)
    {
        if(oldPoint == null)
        {
            StsException.systemError("StsTranslateCoordinates.computeTranslatedPoint() failed." +
                " oldPoints null in computeTranslatedPoints");
            return null;
        }
        if(oldCoor == null)
        {
             StsException.systemError("StsTranslateCoordinates.computeTranslatedPoint() failed." +
                " oldCoor null in computeTranslatedPoints");
            return null;
        }
        if(newCoor == null)
        {
             StsException.systemError("StsTranslateCoordinates.computeTranslatedPoint() failed." +
                " newCoor null in computeTranslatedPoints");
            return null;
        }

        float scaleFactor = newCoor.length/oldCoor.length;

        StsPoint dPoint = StsPoint.subPointsStatic(oldPoint, oldCoor.o);

        float fa = scaleFactor*dPoint.dot(oldCoor.a, 3);
        float fh = scaleFactor*dPoint.dot(oldCoor.h, 3);
        float fv = scaleFactor*dPoint.dot(oldCoor.v, 3);

        float[] v = new float[3];

        for(int i = 0; i < 3; i++)
            v[i] = fa*newCoor.a.v[i] + fh*newCoor.h.v[i] + fv*newCoor.v.v[i] + newCoor.o.v[i];

        return new StsPoint(v);
    }

	// oldCoor includes coordinates of a single old point; returns xyz of point translated to new coordinates
    static public StsPoint computeTranslatedPoint(StsTranslateCoordinates oldCoor, StsTranslateCoordinates newCoor)
    {
        if(oldCoor == null)
        {
             StsException.systemError("StsTranslateCoordinates.computeTranslatedPoint() failed." +
                " oldCoor null in computeTranslatedPoints");
            return null;
        }
        if(newCoor == null)
        {
             StsException.systemError("StsTranslateCoordinates.computeTranslatedPoint() failed." +
                " newCoor null in computeTranslatedPoints");
            return null;
        }

        float scaleFactor = newCoor.length/oldCoor.length;

        float fa = oldCoor.f.v[0];
        float fh = oldCoor.f.v[1];
        float fv = oldCoor.f.v[2];

        float[] v = new float[3];

        for(int i = 0; i < 3; i++)
            v[i] = scaleFactor*(fa*newCoor.a.v[i] + fh*newCoor.h.v[i] + fv*newCoor.v.v[i]) + newCoor.o.v[i];

        return new StsPoint(v);
    }

    public static final void main(String[] args)
    {
		StsTranslateCoordinates oldCoor, newCoor;
		StsPoint newPoint;

		StsPoint oldPoint0 = new StsPoint(0.0f, 0.0f, 0.0f);
		StsPoint oldPoint1 = new StsPoint(1.0f, 1.0f, 1.0f);
		StsPoint oldPoint = new StsPoint(0.5f, 0.5f, 0.5f);
		StsPoint newPoint0 = new StsPoint(1.0f, 1.0f, 1.0f);
		StsPoint newPoint1 = new StsPoint(11.0f, 11.0f, 11.0f);

		oldCoor = new StsTranslateCoordinates(oldPoint0, oldPoint1);
		newCoor = new StsTranslateCoordinates(newPoint0, newPoint1);
		newPoint = StsTranslateCoordinates.computeTranslatedPoint(oldPoint, oldCoor, newCoor);
		newPoint.print();

		oldCoor = new StsTranslateCoordinates(oldPoint0, oldPoint1, oldPoint);
	    newPoint = StsTranslateCoordinates.computeTranslatedPoint(oldCoor, newCoor);
		newPoint.print();
    }
}
