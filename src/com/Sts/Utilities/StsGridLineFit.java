
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Utilities;

import com.Sts.Types.*;

/** Fits a line of the form z = a + bx or z = a + by thru a given set of points
 *  using least-squares. xyDir indicates whether it is a function of x (0) or y (1).
 *  Points may optionally have weights (ignored if null)
 *  a reference point used in computing line sigma for 1 or 2 points (ignored if null).
 */
public class StsGridLineFit
{
    public int nPnts;        /** number of points in line to be fitted       */
    public float x0, y0, s;  /** Coefficients of line equation: y = y0 + s*(x - x0)  */
    public float ySig, sSig; /** Standard deviations for a and b */
    public float lineSig;    /** Standard deviation of the line */
    public float fitQ;       /** Goodness of fit; not currently computed: needs gamma function. */

    static final int ROW = StsParameters.ROW;
    static final int COL = StsParameters.COL;
    static final int PLUS = StsParameters.PLUS;
    static final int MINUS = StsParameters.MINUS;

    /** Least squares line fitting.  Borrowed from
     * "Numerical Recipes in C" 2nd edition, p.665
     */
    public StsGridLineFit(int min, int max, StsPoint2D[] points)
    {
	    int i, n;
	    float wt, t, sxoss;
        float sx, sy, st2;
        float ss, chiSq, chi;
        float sumWtSq;

        nPnts = max - min + 1;
        float[] x = new float[nPnts];
        float[] y = new float[nPnts];
        float[] w = new float[nPnts];

        this.nPnts = nPnts;

        if(nPnts == 1)
        {
            x0 = points[min].x;
            y0 = points[min].y;
            s = 0.0f;
            ySig = 0.0f;
            sSig = 0.0f;
            lineSig = 0.5f*StsParameters.largeFloat;
            fitQ = 0.0f;
            return;
        }

        /** Use the middle point as the reference; set max wt there also */
        int center = (max + min)/2;
        x0 = points[center].x;
        y0 = points[center].y;

        for(i = 0, n = min; i < nPnts; i++, n++)
        {
            x[i] = points[n].x - x0;
            y[i] = points[n].y - y0;
            w[i] = 1.0f/(Math.abs(n - center) + 1);
        }

        if(nPnts == 2)
        {
            s = y[max]/x[max];
            ySig = 0.0f;
            sSig = 0.0f;
            lineSig = 0.5f*Math.abs(y[max]);
            fitQ = 0.0f;
            return;
        }

        sumWtSq = 0.0f;
        sx = 0.0f;
        sy = 0.0f;
        ss = 0.0f;
        for(i = 0; i < nPnts; i++)
        {
            wt = w[i]*w[i];
            sumWtSq += wt;

            ss += wt;
            sx += x[i]*wt;
            sy += y[i]*wt;
        }

        sxoss = sx/ss;

        st2 = 0.0f;
        s = 0.0f;
        for(i = 0; i < nPnts; i++)
        {
            t = (x[i] - sxoss)*w[i];
            st2 += t*t;
            s += t*y[i]*w[i];
        }

        s /= st2;
        y0 = y0 + (sy - sx*s)/ss;

        ySig = (float)Math.sqrt((1.0 + sx*sx/(ss*st2))/ss);
        sSig = (float)Math.sqrt(1.0/st2);

        chiSq = 0.0f;
        for(i = 0, n = min; i < nPnts; i++, n++)
        {
            chi = (points[n].y - y0 - s*(points[n].x - x0))*w[i];
            chiSq += chi*chi;
        }

        lineSig = (float)Math.sqrt(chiSq/sumWtSq);
/*      fitQ = gammq(0.5*(nPnts - 2), 0.5*(chiSq)); */
    }

    public void print()
    {
        System.out.println("y0: " + y0 + " slope: " + s + " ySig: " + ySig +
                           " sSig: " + sSig + " lineSig: " + lineSig);
    }

    public float getSlope()
    {
        return s;
    }

    public float evaluate(int index)
    {
        return y0 + s*((float)index - x0);
    }

    public float getWeight(int index)
    {
        if((float)index == x0)
            return StsParameters.largeFloat;
        else
            return 1.0f/(float)Math.abs((float)index - x0);
    }

    public float getError(int index)
    {
        if((float)index == x0)
            return lineSig;
        else
            return lineSig*Math.abs(index - x0);
    }

    public static void main(String[] args)
    {
        int nPnts = 3;
        StsPoint2D[] points = new StsPoint2D[nPnts];

        points[0] = new StsPoint2D(116.0f, 10439.6f);
        points[1] = new StsPoint2D(115.0f, 10438.0f);
        points[2] = new StsPoint2D(114.0f, 10438.8f);

        StsGridLineFit line = new StsGridLineFit(0, nPnts-1, points);
        line.print();

        float z;
        for(int i = 112; i < 120; i++)
        {
            z = line.evaluate(i);
            System.out.println("i: " + i + " z: " + z);
        }
    }
}
