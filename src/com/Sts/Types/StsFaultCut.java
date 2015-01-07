
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;

/** a set of X-Y points that are used to describe a fault cut on a surface.  It is
    assumed that the points describe a line without any internal segments
    intersecting others within in it.  X-Y values are absolute coordinates.
*/
public class StsFaultCut
{
    private float[][] xy;

    private StsFaultCut beginConnect = null;
    private StsFaultCut endConnect = null;

    public StsFaultCut(double[] x, double[] y)
    {
        int n;

        // error checks
        if (x == null || y == null || x.length != y.length) return;

        StsProject project = StsObject.getCurrentModel().getProject();

        int nValues = x.length;
        xy = new float[nValues][];

        for(n = 0; n < nValues; n++)
            xy[n] = project.getRelativeXY(x[n], y[n]);
    }

    public StsFaultCut(float[][] localXYs)
    {
        this.xy = localXYs;
    }

    public float[][] getXYs() { return xy; }

    public double[][] getAbsoluteXYs()
    {
        if(xy == null || xy.length <= 0) return null;

        int nValues = xy.length;
        double[][] absXYs = new double[nValues][];

        StsProject project = StsObject.getCurrentModel().getProject();

        for(int n = 0; n < nValues; n++)
            absXYs[n] = project.getAbsoluteXYCoordinates(xy[n]);

        return absXYs;
    }

    /** print out */
    public void print()
    {
        System.out.println("StsFaultCut:");
        if (xy == null)
        {
            System.out.println("\tno X-Y values to print.");
            System.out.println();
            return;
        }
        int nPnts = xy.length;
        for (int i=0; i<nPnts; i++)
        {
            System.out.println("\tx[" + i + "] = " + xy[i][0] + "\ty[" + i + "] = "
                    + xy[i][1]);
        }
        System.out.println();
    }

    public boolean clipToGrid(StsModelSurface surface)
    {
        int nOK = 0;
        int size = xy.length;
        float[][] okXY = new float[size][2];
        float xMin = surface.getXMin();
        float xMax = surface.getXMax();
        float yMin = surface.getYMin();
        float yMax = surface.getYMax();

        for( int i=0; i<size; i++ )
        {
            if( xy[i][0] >= xMin &&
                xy[i][0] <= xMax &&
                xy[i][1] >= yMin &&
                xy[i][1] <= yMax )
            {
                okXY[nOK++] = xy[i];
            }
        }
        if(nOK == size)
        {
            okXY = null;
            return true;
        }
        if(nOK == 0)
        {
            okXY = null;
            return false;
        }
        else
        {
            System.arraycopy(okXY, 0, xy, 0, nOK);
            return true;
        }
    }

    public void testPrint(int k)
    {
        if (xy == null) return;
        System.out.println();
        float z = StsParameters.nullValue;
        int nPnts = xy.length;
        System.out.println(xy[0][0] + "\t" + xy[0][1] + "\t" + z + "\t6\t" + k);
        for (int i=1; i<nPnts-1; i++)
        {
            System.out.println(xy[i][0] + "\t" + xy[i][1] + "\t" + z + "\t7\t" + k);
        }
        if (nPnts>1) System.out.println(xy[nPnts-1][0] + "\t" + xy[nPnts-1][1] + "\t" + z + "\t8\t" + k);

        System.out.println();
    }
}
