
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Utilities;

import com.Sts.Types.*;

public class StsBezier
{
//	static Matrix4d bezierMatrix = new Matrix4d( -1.0,  3.0, -3.0,  1.0,  3.0, -6.0,  3.0,  0.0,  -3.0,  3.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0);
//	static Matrix4d precisionMatrix = null;

//	Matrix4d fdMatrix; // forward difference matrix for fast interpolation

    static public StsPoint evalBezierCurve(StsPoint[] controlPoints, int degree, float u)
	{
    	int i, j;

        StsPoint[] points  = new StsPoint[degree+1];

        for(i = 0; i <= degree; i++)
            points[i] = new StsPoint(controlPoints[i]);

        for(i = 1; i <= degree; i++)
            for(j = 0; j <= degree-i; j++)
                points[j].interpolatePoints(points[j], points[j+1], u);

        return points[0];
    }

    static public StsPoint[] computeXYZorTLine(StsPoint[] controlPoints, int pointsPerVertex, boolean isDepth)
    {
        try
        {
            if(controlPoints == null || controlPoints.length < 2) return controlPoints;

            if(!isDepth)
            {
                for(int n = 0; n < controlPoints.length; n++)
                {
                    float[] xyzmt = controlPoints[n].v;
                    float tempZ = xyzmt[2];
                    xyzmt[2] = xyzmt[4];
                    xyzmt[4] = tempZ;
                }
            }

            StsPoint[] slopes = computeXYZSlopes(controlPoints);
            if(slopes == null) return null;
            if(!isDepth)
            {
                for(int n = 0; n < controlPoints.length; n++)
                {
                    float[] xyzmt = controlPoints[n].v;
                    float tempT = xyzmt[2];
                    xyzmt[2] = xyzmt[4];
                    xyzmt[4] = tempT;
                }
            }
            return computeXYZPoints(controlPoints, slopes, pointsPerVertex);
        }
        catch( Exception e )
        {
            StsException.outputException("System Error in StsBezier.computeXYZLine()", e, StsException.WARNING);
        }

        return null;
    }

    static public StsPoint[] computeXYZLine(StsPoint[] controlPoints, int pointsPerVertex)
    {
        try
        {
			if(controlPoints == null || controlPoints.length < 2) return controlPoints;

        	StsPoint[] slopes = computeXYZSlopes(controlPoints);
            if(slopes == null)
                return null;
            else
                return computeXYZPoints(controlPoints, slopes, pointsPerVertex);
        }
        catch( Exception e )
        {
            StsException.outputException("System Error in StsBezier.computeXYZLine()", e, StsException.WARNING);
        }

        return null;
    }

  	static public StsPoint[] computeXYZSlopes(StsPoint[] controlPoints) throws StsException
    {
    	float len1, len2, f;

        if(controlPoints == null) return null;
        int nControlPoints = controlPoints.length;
        if(nControlPoints <= 1) return null;

        StsPoint[] segmentSlopes = new StsPoint[nControlPoints-1];
        StsPoint[] slopes = new StsPoint[nControlPoints];

    	for(int i = 1; i < nControlPoints; i++)
        {
       		StsPoint segmentSlope = new StsPoint(4);
			segmentSlope.subPoints(controlPoints[i], controlPoints[i-1]);
            segmentSlope.v[3] = segmentSlope.normalizeReturnLength();
   			segmentSlopes[i-1] = segmentSlope;
        }

        /* Slope at end-points is defined by slope with nearest point */
        /* intermediate points have weighted slopes					  */

        slopes[0] = segmentSlopes[0];

		len2 = segmentSlopes[0].v[3];
		for(int i = 1; i < nControlPoints-1; i++)
		{
			len1 = len2;
			len2 = segmentSlopes[i].v[3];
			f = len1/(len1 + len2);
            slopes[i] = new StsPoint(4);
			slopes[i].interpolatePoints(3, segmentSlopes[i-1], segmentSlopes[i], f);
			slopes[i].v[3] = segmentSlopes[i].v[3];
		}
        slopes[nControlPoints-1] = segmentSlopes[nControlPoints-2];

        segmentSlopes = null;

		return slopes;
    }

    static public StsPoint[] computeXYZPoints(StsPoint[] controlPoints, StsPoint[] slopes, int pointsPerVertex) throws StsException
    {
        StsPoint[] points;
        float length;
        int degree = 3;
        int dim = 3;
        StsPoint[] cp = new StsPoint[degree+1];
        float u, du;
        int i, n;

        int nControlPoints = controlPoints.length;
        if(nControlPoints < 1) return null;

        if(nControlPoints == 1 || pointsPerVertex == 1)
            return controlPoints;

        int nPoints = (nControlPoints-1)*pointsPerVertex + 1;
        points = new StsPoint[nPoints];

        int nPnt = 0;

        /** Insert vertexPoint into points and then pointsPerVertex-1 interpolated points */

    	for(i = 0; i < nControlPoints-1; i++)
        {
            cp[0] = controlPoints[i];
            cp[3] = controlPoints[i+1];

            length = slopes[i].v[3]/3.0f;
            cp[1] = StsPoint.multByConstantAddPointStatic(slopes[i], length, cp[0]);
            cp[2] = StsPoint.multByConstantAddPointStatic(slopes[i+1], -length, cp[3]);

            du = 1.0f/pointsPerVertex;
            u = 0.0f;

            points[nPnt++] = cp[0];

            for(n = 1; n < pointsPerVertex; n++)
            {
                u += du;
                points[nPnt++] = StsBezier.evalBezierCurve(cp, degree, u);
            }
        }

        /** Insert last point into points */
        points[nPnt] = cp[3];

        return points;
	}

    public static float[][] createCubicControlPoints1D(float value0, float value1, float slope0, float slope1)
    {
        float[][] controlPoints = new float[4][2];
        controlPoints[0] = new float[] { value0, 0.0f };
        controlPoints[1] = new float[] { value0 + slope0/3.0f, 0.33333f };
        controlPoints[2] = new float[] { value1 - slope1/3.0f, 0.66666f };
        controlPoints[3] = new float[] { value1, 1.0f };
        return controlPoints;
    }

    static public float[][] subdivide(int degree, float[] in, float t)
    {
        float t1 = 1.0f - t;

        float[] outLeft = new float[degree+1];
        float[] outRite = new float[degree+1];

        for(int i = 0; i <= degree; i++)
            outRite[i] = in[i];

        for(int r = 1; r <= degree; r++)
            for(int i = 0; i <= degree - r; i++)
                outRite[i] = t1*outRite[i] + t*outRite[i+1];

        t = t1;
        t1 = 1.0f - t;

        for(int i = 0; i <= degree; i++)
            outLeft[degree - i] = in[i];

        for(int r = 1; r <= degree; r++)
            for(int i = 0; i <= degree - r; i++)
                outLeft[i] = t1*outLeft[i] + t*outLeft[i+1];

        return new float[][] { outLeft, outRite };
    }

	/** requires that measured length along path is in the measured depth field  (4th field). */
	static public StsPoint[] computeXYZLineSlopes(StsPoint[] points)
	{
		StsPoint slope1, slope2;
		float len1, len2, f;

		int nPoints = points.length;
		if (nPoints <= 1)
		{
			return null;
		}
		int nCoors = points[0].v.length;
		StsPoint[] slopes = new StsPoint[nPoints];
		for (int i = 0; i < nPoints - 1; i++)
		{
			slopes[i] = new StsPoint(nCoors);
			slopes[i].subPoints(points[i+1], points[i]);
		}
		slopes[nPoints-1] = new StsPoint(slopes[nPoints-2]);

		/* Slope at end-points is defined by slope with nearest point */
		/* intermediate points have weighted slopes					  */

//        vertexSlopes[0].print("vertexSlope: ", 0);

		len2 = slopes[0].getM();
		slope2 = slopes[0];

		for (int i = 1; i < nPoints - 1; i++)
		{
			len1 = len2;
			slope1 = slope2;

			len2 = slopes[i].getM();
			slope2 = slopes[i];

			f = len1 / (len1 + len2);
			slopes[i].interpolatePoints(slope1, slope2, f);
//         	//vertexSlopes[i].print("vertexSlope: ", i);
		}
//		slopes[nPoints - 1] = new StsPoint(slopes[nPoints - 2]);
		for(int i = 0; i < nPoints; i++)
		{
			normalizeByIntervalLength(slopes[i]);
		}
//      //vertexSlopes[nVertices-1].print("vertexSlope: ", nVertices-1);
		return slopes;
	}

	/** We expect the distance from this point to the next to be in the measured depth position.
	 *  Normalize the other values by this distance value.
	 */
	static private void normalizeByIntervalLength(StsPoint point)
	{
		float length = point.getM();
		if(length <= 0.0f) return;
		float[] v = point.v;
		for(int n = 0; n < v.length; n++)
			v[n] /= length;
    }

	/** Given a series of XYZM points and slopes where M is the parameter, Bezier cubic interpolate
	 *  over the parameter range from pMin to pMax where pMax is defined by pInc and nNewPoints */
	static public StsPoint[] computeXYZPoints(StsPoint[] points, StsPoint[] slopes, float pMin, float pInc, int nNewPoints)
	{
		float u;

		int nPoints = points.length;
		if (nPoints < 1)return null;

		if (nPoints == 1 || slopes == null)return null;

		StsPoint[] newPoints = new StsPoint[nNewPoints];
		float p = pMin;

		StsPoint point0 = points[0];
		StsPoint point1 = points[1];
		StsPoint slope0 = slopes[0];
		StsPoint slope1 = slopes[1];
		StsPoint[] controlPoints = getControlPoints(point0, point1, slope0, slope1);
		float p0 = point0.getM();
		float p1 = point1.getM();
		int i = 0;
		for (int n = 0; n < nNewPoints; n++, p += pInc)
		{
			if (p < p0)
			{
				u = p - p0;
				newPoints[n] = new StsPoint(4);
				newPoints[n].multByConstantAddPoint(slope0, u, point0);
				continue;
			}
			else if (p < p1)
			{
				u = (p - p0) / (p1 - p0);
				newPoints[n] = StsBezier.evalBezierCurve(controlPoints, 3, u);
				continue;
			}
			while (p >= p1 && i < nPoints - 1)
			{
				point0 = point1;
				slope0 = slope1;
				i++;
				point1 = points[i];
				slope1 = slopes[i];
				p0 = p1;
				p1 = point1.getM();
				controlPoints = getControlPoints(point0, point1, slope0, slope1);
			}
			if (p < p1)
			{
				u = (p - p0) / (p1 - p0);
				newPoints[n] = StsBezier.evalBezierCurve(controlPoints, 3, u);
			}
			else
			{
				u = p - p1;
				newPoints[n] = new StsPoint(4);
				newPoints[n].multByConstantAddPoint(slope1, u, point1);
			}
		}
		return newPoints;
	}

	static private StsPoint[] getControlPoints(StsPoint point0, StsPoint point1, StsPoint slope0, StsPoint slope1)
	{
		StsPoint[] cp = new StsPoint[4];
		cp[0] = point0;
		cp[3] = point1;
		float len = cp[3].getM() - cp[0].getM();

		float f = len /3.0f;
		cp[1] = new StsPoint(4);
		cp[1].multByConstantAddPoint(slope0, f, point0);

		f = -len /3.0f;
		cp[2] = new StsPoint(4);
		cp[2].multByConstantAddPoint(slope1, f, point1);
		return cp;
	}
/*
	public void makeCurveFdMatrix(StsPoint[] controlPoints, int nCurvePoints)
	{
		fdMatrix = new Matrix4d();
		for(int row = 0; row < 4; row++)
		{
			for(int col = 0; col < 3; col++)
				fdMatrix.setElement(row, col, controlPoints[row].v[col]);
			fdMatrix.setElement(row, 3, 1.0);
		}

		if(precisionMatrix == null) makePrecisionMatrix(1.0/(nCurvePoints-1));
		Matrix4d matrix = new Matrix4d(precisionMatrix);
		matrix.mul(bezierMatrix);
		matrix.mul(fdMatrix);
		fdMatrix = matrix;
	}

	static void makePrecisionMatrix(double rn)
	{
		double rn2 = rn*rn;
		double rn3 = rn*rn2;

		precisionMatrix = new Matrix4d();
		precisionMatrix.setIdentity();

		precisionMatrix.setElement(0, 0, 6.0*rn3);
		precisionMatrix.setElement(1, 0, 6.0*rn3);
		precisionMatrix.setElement(2, 0, rn3);

		precisionMatrix.setElement(1, 1, 2.0*rn2);
		precisionMatrix.setElement(2, 1, rn2);
		precisionMatrix.setElement(2, 2, rn);
	}

	public StsPoint[] generateBezierCurve(StsPoint[] controlPoints, int nCurvePoints)
	{
		StsPoint[] curvePoints = new StsPoint[nCurvePoints];
		double rscale;
		int i, j, n;

		makeCurveFdMatrix(controlPoints, nCurvePoints);
		for(n = 0; n < nCurvePoints; n++)
		{
			curvePoints[n] = new StsPoint(3);
			rscale  = 1.0/fdMatrix.getElement(3, 3);

			curvePoints[n].v[0] = (float)(fdMatrix.getElement(3, 0)*rscale);
			curvePoints[n].v[1] = (float)(fdMatrix.getElement(3, 1)*rscale);
			curvePoints[n].v[2] = (float)(fdMatrix.getElement(3, 2)*rscale);

			for(i = 3; i > 0; i--)
			{
				for(j = 0; j <= 3; j++)
				{
					double value = fdMatrix.getElement(i, j);
					double dValue = fdMatrix.getElement(i - 1, j);
					fdMatrix.setElement(i, j, value + dValue);
				}
			}
		}
		return curvePoints;
	}

	static public void main(String[] args)
	{
		StsBezier bezier = new StsBezier();
		StsPoint[] controlPoints = new StsPoint[4];
		controlPoints[0] = new StsPoint(0.0, 0.0, 0.0);
		controlPoints[1] = new StsPoint(1.0, 1.0, 0.0);
		controlPoints[2] = new StsPoint(2.0, 1.0, 0.0);
		controlPoints[3] = new StsPoint(3.0, 0.0, 0.0);
		int nCurvePoints = 1001;
		StsPoint[] curvePoints = bezier.generateBezierCurve(controlPoints, nCurvePoints);
		float du = 1.0f/(nCurvePoints-1);
		float u = 0.0f;
		for(int n = 0; n < nCurvePoints; n++, u += du)
		{
			StsPoint evalPoint = evalBezierCurve(controlPoints, 3, u);
			curvePoints[n].print("point     ", n);
			evalPoint.print("  compare ", n);
		}
	}
*/
}
