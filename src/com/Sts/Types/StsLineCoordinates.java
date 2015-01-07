
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.Utilities.*;
/** Defines a coordinate system for a line between two points */

public class StsLineCoordinates
{
    /**  	@param a is unit normal along line between first and last points
      *		@param h is unit normal orthogonal to a and horizontal thru first point
      *		@param v is unit normal orthogonal to a and h.
      *     @param o is origin (first point).
	  *     @param f is point in a, h, v coordinate system. 4th argument is control parameter for interpolation.
      *     @param length is distance between origin and last point
      */

    public StsPoint a, h, v, o;
	public StsPoint f;
    public float length;

    public StsLineCoordinates(StsPoint firstPoint, StsPoint lastPoint)
    {
        StsPoint firstXYZPoint = firstPoint;
        StsPoint lastXYZPoint = lastPoint;
        initialize(firstXYZPoint, lastXYZPoint);
    }

    public StsLineCoordinates(StsPoint firstPoint, StsPoint lastPoint, boolean isDepth)
    {
        StsPoint firstXYZPoint = firstPoint.getXYZorTPoint(isDepth);
        StsPoint lastXYZPoint = lastPoint.getXYZorTPoint(isDepth);
        initialize(firstXYZPoint, lastXYZPoint);
    }

    public StsLineCoordinates(StsPoint firstPoint, StsPoint lastPoint, StsPoint point, float z)
    {
        initialize(firstPoint, lastPoint);
		StsPoint dPoint = StsPoint.subPointsStatic(point, o, 3);
		float pointLength = dPoint.normalizeReturnLength();
		float scaleFactor = pointLength/length;
		f = new StsPoint(4);
        f.v[0] = scaleFactor*dPoint.dot(a, 3);
        f.v[1] = scaleFactor*dPoint.dot(h, 3);
        f.v[2] = scaleFactor*dPoint.dot(v, 3);
		f.v[3] = z;
    }

    private void initialize(StsPoint firstXYZPoint, StsPoint lastXYZPoint)
    {
        o = firstXYZPoint;
        a = StsPoint.subPointsStatic(lastXYZPoint, firstXYZPoint);

//        o = new StsPoint(firstPoint);
//        a = new StsPoint();
//        a.subPoints(lastPoint, firstPoint);
        length = a.normalizeReturnLength();
        h = new StsPoint(a.v[1],  -a.v[0], 0.0f);
        h.normalizeXYZ();

        v = new StsPoint(3);
        v.crossProduct(a, h);
    }

    /** With a point defined in this lineCoordinate system, compute location of actual point. */
	public StsPoint computeTranslatedPoint(StsPoint point)
	{
		float fa = point.v[0];
		float fh = point.v[1];
		float fv = point.v[2];
		StsPoint newPoint = new StsPoint(4);
		for(int i = 0; i < 3; i++)
			// newPoint.v[i] = fa*a.v[i] + fh*h.v[i] + fv*v.v[i] + o.v[i];
            newPoint.v[i] = length*(fa*a.v[i] + fh*h.v[i] + fv*v.v[i]) + o.v[i];
        return newPoint;
	}

    static public void computeTranslatedPoints(StsPoint[] oldPoints, StsLineCoordinates oldCoor,
                                               StsPoint[] newPoints, StsLineCoordinates newCoor)
        throws StsException
    {
        StsPoint dPoint = new StsPoint(3);

        if(oldPoints == null)
        {
            StsException.systemError("StsLineCoordinates.computeTranslatedPoint() failed." +
                " oldPoints null in computeTranslatedPoints");
            return;
        }
        if(oldCoor == null)
        {
            StsException.systemError("StsLineCoordinates.computeTranslatedPoint() failed." +
                " oldCoor null in computeTranslatedPoints");
            return;
        }
        if(newPoints == null)
        {
            StsException.systemError("StsLineCoordinates.computeTranslatedPoint() failed." +
                " newPoints null in computeTranslatedPoints");
            return;
        }
        if(newCoor == null)
        {
            StsException.systemError("StsLineCoordinates.computeTranslatedPoint() failed." +
                " newCoor null in computeTranslatedPoints");
            return;
        }

        float scaleFactor = newCoor.length/oldCoor.length;

        /** now compute new internal points by scaling from old points */
        for(int n = 0; n < oldPoints.length; n++)
        {
            dPoint.subPoints(oldPoints[n], oldCoor.o);

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

    static public void computeTranslatedPoints(StsLineCoordinates[] oldCoors,
                                               StsLineCoordinates newCoor, StsPoint[] newPoints)
        throws StsException
    {
        if(oldCoors == null)
        {
            StsException.systemError("StsLineCoordinates.computeTranslatedPoint() failed." +
                " oldCoors null in computeTranslatedPoints");
            return;
        }
        if(newPoints == null)
        {
            StsException.systemError("StsLineCoordinates.computeTranslatedPoint() failed." +
                " newPoints null in computeTranslatedPoints");
            return;
        }
        if(newCoor == null)
        {
            StsException.systemError("StsLineCoordinates.computeTranslatedPoint() failed." +
                " newCoor null in computeTranslatedPoints");
            return;
        }

		int nPoints = oldCoors.length;
        /** now compute new internal points by scaling from old points */
        for(int n = 0; n < nPoints; n++)
        {
            float scaleFactor = newCoor.length/oldCoors[n].length;
		    float[] oldV = oldCoors[n].f.v;
			float fa = oldV[0];
			float fh = oldV[1];
			float fv = oldV[2];

			float[] newA = newCoor.a.v;
			float[] newH = newCoor.h.v;
			float[] newV = newCoor.v.v;
			float[] newO = newCoor.o.v;
			float[] v = new float[4];

			for(int i = 0; i < 3; i++)
				v[i] = scaleFactor*(fa*newA[i] + fh*newH[i] + fv*newV[i]) + newO[i];

            if(newPoints[n] == null)
                newPoints[n] = new StsPoint(v);
            else
                newPoints[n].setValues(v);
        }
    }

    static public StsPoint[] computeTranslatedPoints(StsPoint[] oldPoints, StsLineCoordinates oldCoor,
                                                     StsPoint firstPoint, StsPoint lastPoint)
        throws StsException
    {
        StsLineCoordinates newCoor = new StsLineCoordinates(firstPoint, lastPoint);

        StsPoint[] newPoints = new StsPoint[oldPoints.length];

        computeTranslatedPoints(oldPoints, oldCoor, newPoints, newCoor);

        return newPoints;
    }

	// computes oldCoordinates of oldPoint and returns xyz of point translated to newCoorinates
    static public StsPoint computeTranslatedPoint(StsPoint oldPoint, StsLineCoordinates oldCoor,
                                                  StsLineCoordinates newCoor)
    {
        if(oldPoint == null)
        {
            StsException.systemError("StsLineCoordinates.computeTranslatedPoint() failed." +
                " oldPoints null in computeTranslatedPoints");
            return null;
        }
        if(oldCoor == null)
        {
             StsException.systemError("StsLineCoordinates.computeTranslatedPoint() failed." +
                " oldCoor null in computeTranslatedPoints");
            return null;
        }
        if(newCoor == null)
        {
             StsException.systemError("StsLineCoordinates.computeTranslatedPoint() failed." +
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
    static public StsPoint computeTranslatedPoint(StsLineCoordinates oldCoor, StsLineCoordinates newCoor)
    {
        if(oldCoor == null)
        {
             StsException.systemError("StsLineCoordinates.computeTranslatedPoint() failed." +
                " oldCoor null in computeTranslatedPoints");
            return null;
        }
        if(newCoor == null)
        {
             StsException.systemError("StsLineCoordinates.computeTranslatedPoint() failed." +
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
		StsPoint[] firstColVertexPoints = new StsPoint[3];
		firstColVertexPoints[0] = new StsPoint(0.0f, 0.0f, 0.0f, 0.0f);
		firstColVertexPoints[1] = new StsPoint(0.0f, 0.0f, 1.0f, 1.0f);
		firstColVertexPoints[2] = new StsPoint(0.0f, 0.0f, 2.0f, 2.0f);

		StsPoint[] lastColVertexPoints = new StsPoint[3];
	    lastColVertexPoints[0] = new StsPoint(1.0f, 0.0f, 0.0f, 0.0f);
	    lastColVertexPoints[1] = new StsPoint(1.0f, 1.0f, 1.0f, 1.0f);
	    lastColVertexPoints[2] = new StsPoint(0.0f, 1.0f, 2.0f, 2.0f);

	    StsLineCoordinates[] lineCoordinates = new StsLineCoordinates[3];
		for(int n = 0; n < 3; n++)
			lineCoordinates[n] = new StsLineCoordinates(firstColVertexPoints[n], lastColVertexPoints[n]);

	    int nRows = 21;
		float zMin = 0.0f;
		float zInc = 0.1f;

	    StsPoint[] firstSlopes = StsBezier.computeXYZLineSlopes(firstColVertexPoints);
	    StsPoint[] firstColPoints = StsBezier.computeXYZPoints(firstColVertexPoints, firstSlopes, zMin, zInc, nRows);

	    StsPoint[] lastSlopes = StsBezier.computeXYZLineSlopes(lastColVertexPoints);
	    StsPoint[] lastColPoints = StsBezier.computeXYZPoints(lastColVertexPoints, lastSlopes, zMin, zInc, nRows);

	    StsPoint[] lineCoordinatePoints = new StsPoint[3];
		lineCoordinatePoints[0] = new StsPoint(0.5f, 0.0f, 0.0f, 0.0f);
		lineCoordinatePoints[1] = new StsPoint(0.5f, 0.0f, 0.0f, 1.0f);
		lineCoordinatePoints[2] = new StsPoint(0.5f, 0.0f, 0.0f, 2.0f);

	    StsPoint[] lineCoordinateSlopes = StsBezier.computeXYZLineSlopes(lineCoordinatePoints);
		StsPoint[] lineCoordinateColPoints = StsBezier.computeXYZPoints(lineCoordinatePoints, lineCoordinateSlopes, zMin, zInc, nRows);

		StsPoint[] colPoints = new StsPoint[nRows];
		for (int row = 0; row < nRows; row++)
		{
			StsLineCoordinates rowCoordinates = new StsLineCoordinates(firstColPoints[row], lastColPoints[row]);
			colPoints[row] = rowCoordinates.computeTranslatedPoint(lineCoordinateColPoints[row]);
		}
    }
}
