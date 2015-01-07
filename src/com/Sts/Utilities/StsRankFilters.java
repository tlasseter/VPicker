//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Utilities;

import com.Sts.Actions.Wizards.VirtualVolume.*;


public class StsRankFilters
{
    public static final String MEAN = StsFilterVVolumePanel.MEAN;
    public static final String MEDIAN = StsFilterVVolumePanel.MEDIAN;
    public static final String MAX = StsFilterVVolumePanel.MAX;
    public static final String MIN = StsFilterVVolumePanel.MIN;
    public static final String VARIANCE = StsFilterVVolumePanel.VARIANCE;

    public static final int XDIR = StsParameters.XDIR;
    public static final int YDIR = StsParameters.YDIR;
    public static final int ZDIR = StsParameters.ZDIR;

    static public float[] rankFloat3D(int dir, int nPlaneRows, int nPlaneCols, int nPlaneMin, int nPlaneMax,
                                      float[][] inPlaneData, String filterName, int xSize, int ySize, int zSize)
    {
        int nPlanes = nPlaneMax - nPlaneMin + 1;
        int rowSize = getRowSize(dir, xSize, zSize);
        int colSize = getColSize(dir, xSize, ySize);

        if( filterName.equals(MEDIAN))
            return applyMedianFilter(inPlaneData, nPlanes, nPlaneRows, nPlaneCols, rowSize, colSize);
        else if( filterName.equals(MEAN))
            return applyMeanFilter(inPlaneData, nPlanes, nPlaneRows, nPlaneCols, rowSize, colSize);
        else if( filterName.equals(MIN))
            return applyMinFilter(inPlaneData, nPlanes, nPlaneRows, nPlaneCols, rowSize, colSize);
        else if( filterName.equals(MAX))
            return applyMaxFilter(inPlaneData, nPlanes, nPlaneRows, nPlaneCols, rowSize, colSize);
        else if( filterName.equals(VARIANCE))
            return applyVarianceFilter(inPlaneData, nPlanes, nPlaneRows, nPlaneCols, rowSize, colSize);

        return null;
    }

    static private int getRowSize( int dir, int xSize, int zSize)
    {
        switch(dir)
        {
            case XDIR:
                return zSize;
            case YDIR:
                return zSize;
            case ZDIR:
                return xSize;
            default:
                return 0;
        }
    }

    static private int getColSize( int dir, int xSize, int ySize)
    {
        switch(dir)
        {
            case XDIR:
                return xSize;
            case YDIR:
                return ySize;
            case ZDIR:
                return ySize;
            default:
                return 0;
        }
    }

    static private float[] getVolumeValues( float[][] inPlaneData, int row, int col, int nPlanes, int nPlaneRows,
                                            int nPlaneCols, int rowSize, int colSize)
    {
        int u1 = Math.max(col - colSize/2, 0);
        int u2 = Math.min(col + colSize/2, nPlaneCols - 1);
        int v1 = Math.max(row - rowSize/2, 0);
        int v2 = Math.min(row + rowSize/2, nPlaneRows - 1);

        int nElements = (u2 - u1 + 1)*(v2 - v1 + 1);
        float [] volumeValues = new float[nPlanes*nElements];

        int count = 0;
        for( int plane=0; plane<nPlanes; plane++)
        {
            for( int v=v1; v<=v2; v++)
            {
                for( int u=u1; u<=u2; u++)
                    volumeValues[count++] = inPlaneData[plane][u + v*nPlaneCols];
            }
        }
        return volumeValues;
    }

    static private float[] applyMedianFilter( float[][]inPlaneData, int nPlanes, int nPlaneRows, int nPlaneCols,
                                              int rowSize, int colSize)
    {
       float [] outPlaneData = new float[nPlaneRows*nPlaneCols];
       float [] values;

       for( int row=0; row<nPlaneRows; row++)
        {
            for( int col=0; col<nPlaneCols; col++)
            {
                values = getVolumeValues( inPlaneData, row, col, nPlanes, nPlaneRows, nPlaneCols, rowSize, colSize);
                outPlaneData[col+row*nPlaneCols] = findMedian(values);
            }
        }
        return outPlaneData;
    }

    static private float findMedian(float[] a)
    {
        final int nValues = a.length;
        final int nv1b2 = (nValues-1)/2;
        int i,j;
        int l=0;
        int m=nValues-1;
        float med=a[nv1b2];
        float dum ;

        while (l<m)
        {
            i=l ;
            j=m ;
            do {
                while (a[i]<med) i++ ;
                while (med<a[j]) j-- ;
                dum=a[j];
                a[j]=a[i];
                a[i]=dum;
                i++ ; j-- ;
            } while ((j>=nv1b2) && (i<=nv1b2)) ;
            if (j<nv1b2) l=i ;
            if (nv1b2<i) m=j ;
            med=a[nv1b2] ;
        }
        return med ;
    }

    static private float[] applyMinFilter( float[][]inPlaneData, int nPlanes, int nPlaneRows, int nPlaneCols,
                                           int rowSize, int colSize)
    {
       float [] outPlaneData = new float[nPlaneRows*nPlaneCols];
       float [] values;

       for( int row=0; row<nPlaneRows; row++)
        {
            for( int col=0; col<nPlaneCols; col++)
            {
                values = getVolumeValues( inPlaneData, row, col, nPlanes, nPlaneRows, nPlaneCols, rowSize, colSize);
                outPlaneData[col+row*nPlaneCols] = findMin(values);
            }
        }
        return outPlaneData;
    }

    static private float findMin(float[] values)
    {
        float min = values[0];
        for (int i=1; i<values.length; i++)
            if (values[i]<min)
                min = values[i];
        return min;
    }

    static private float[] applyMaxFilter( float[][]inPlaneData, int nPlanes, int nPlaneRows, int nPlaneCols,
                                           int rowSize, int colSize)
    {
       float [] outPlaneData = new float[nPlaneRows*nPlaneCols];
       float [] values;

       for( int row=0; row<nPlaneRows; row++)
        {
            for( int col=0; col<nPlaneCols; col++)
            {
                values = getVolumeValues( inPlaneData, row, col, nPlanes, nPlaneRows, nPlaneCols, rowSize, colSize);
                outPlaneData[col+row*nPlaneCols] = findMax(values);
            }
        }
        return outPlaneData;
    }

    static private float findMax(float[] values)
    {
        float max = values[0];
        for (int i=1; i<values.length; i++)
            if (values[i]>max)
                max = values[i];
        return max;
    }

    static private float[] applyMeanFilter( float[][]inPlaneData, int nPlanes, int nPlaneRows, int nPlaneCols,
                                            int rowSize, int colSize)
    {
       float [] outPlaneData = new float[nPlaneRows*nPlaneCols];
       float [] values;

       for( int row=0; row<nPlaneRows; row++)
        {
            for( int col=0; col<nPlaneCols; col++)
            {
                values = getVolumeValues( inPlaneData, row, col, nPlanes, nPlaneRows, nPlaneCols, rowSize, colSize);
                outPlaneData[col+row*nPlaneCols] = findMean(values);
            }
        }
        return outPlaneData;
    }

    static private float findMean(float[] values)
    {
        float sum = values[0];
        for (int i=1; i<values.length; i++)
            sum += values[i];
        return (sum/values.length);
    }

    static private float[] applyVarianceFilter( float[][]inPlaneData, int nPlanes, int nPlaneRows, int nPlaneCols,
                                                int rowSize, int colSize)
    {
       float [] outPlaneData = new float[nPlaneRows*nPlaneCols];
       float [] values;

       for( int row=0; row<nPlaneRows; row++)
        {
            for( int col=0; col<nPlaneCols; col++)
            {
                values = getVolumeValues( inPlaneData, row,col, nPlanes, nPlaneRows, nPlaneCols, rowSize, colSize);
                outPlaneData[col+row*nPlaneCols] = findVariance(values);
            }
        }
        return outPlaneData;
    }

    static private float findVariance(float[] values)
    {
        double v, sum=0.0, sum2=0.0;
        float min = findMin(values);
        int n = values.length;
        for (int i=1; i<n; i++)
        {
            v = values[i] - min;
            sum += v;
            sum2 += v*v;
        }
        double variance = (n*sum2-sum*sum)/n;
        return (float)variance;
    }
}
