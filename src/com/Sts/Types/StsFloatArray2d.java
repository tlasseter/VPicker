
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

// Allocates a single-dimensioned array of points where points have
// a dimension rank: array[nPnts][rank].

// These arrays can be either fixed or growable.

// If fixed, call them with incSize of 0.  size is set to initSize
// and array is flagged as fixed. Points are set using set<type>
// where <type> is the class of the point passed in.

// If growable, call them with positive incSize.  size is set to 0
// and points are added to the array with the add method.

package com.Sts.Types;

import com.Sts.Utilities.*;

public class StsFloatArray2d
{
	private boolean growable;		// this is growable array (not fixed)
    public int nPntsUsed;			// no of points currently used
	public int nPntsAlloc;			// number of allocated points
    public int dimPnt;              // dimension of points
    private int incPnts;            // 0 if fixed; number to increment by if growable
    public float array[][] = null;  // allocated 2-dimensional array

	public StsFloatArray2d(int nPnts, int incPnts, int dimPnt) throws StsException
	{
    	growable = (incPnts > 0);

        if(growable)
        	nPntsUsed = 0;
        else
        	nPntsUsed = nPnts;

    	this.incPnts = incPnts;
        this.dimPnt = dimPnt;

        try
        {
        	array = new float[nPnts][dimPnt];
            nPntsAlloc = nPnts;
        }
        catch (Exception e)
    	{
    		throw new StsException("StsFloatArray2d", "Error allocating 2d float array of size : ", nPnts, " by ", dimPnt );
        }
	}

	public StsFloatArray2d(int nPnts, int dimPnt) throws StsException
	{
    	// This constructTraceAnalyzer makes a fixed array

    	growable = false;
       	nPntsUsed = nPnts;
        incPnts = 0;
        nPntsUsed = nPnts;
        this.dimPnt = dimPnt;

        try
        {
        	array = new float[nPnts][dimPnt];
        	nPntsAlloc = nPnts;
        }
        catch (Exception e)
    	{
    		throw new StsException("StsFloatArray2d.constructor1", "Can't allocate 2d float array of size : ", nPnts, " by ", dimPnt );
        }
	}

    /** accessors */
    public void setValues(float[][] values) { array = values; }
    public float[][] getValues() { return array; }

    public float[] elementAt(int index) throws StsException
    {
        if(index >= nPntsUsed)
        	throw new StsException("StsFloatArray2d.constructor2", "Requested index: ", index, ">= no current points: ", nPntsUsed);
        else
        	return array[index];
    }

    public void add(float[] points) throws StsException
    {
    	if(!growable)
        	throw new StsException("StsFloatArray2d.add(float[])", "Can't use add method for fixed array.");

        if(nPntsAlloc == nPntsUsed) resize(nPntsAlloc+incPnts);

        if(points.length > dimPnt)
        {
        	throw new StsException("StsFloatArray2d.add(float[])", "Point dimension: ",
            					   points.length, "> array dimension: ", dimPnt);
        }
        else
        {
        	for(int n = 0; n < dimPnt; n++)
        		array[nPntsUsed][n] = points[n];

        	nPntsUsed++;
        }
    }

    public void add(StsPoint point) throws StsException
    {
    	if(!growable)
        	throw new StsException("StsFloatArray2d.add(StsPoint)", "Can't use add method for fixed array.");

        if(nPntsAlloc == nPntsUsed) resize(nPntsAlloc+incPnts);

        if(point.v.length > dimPnt)
        {
        	throw new StsException("StsFloatArray2d.add(StsPoint)", "Point dimension: ",
            					   point.v.length, "> array dimension: ", dimPnt);
        }
        else
        {
        	for(int n = 0; n < dimPnt; n++)
        		array[nPntsUsed][n] = point.v[n];

        	nPntsUsed++;
        }
    }

    public void add(float v0, float v1, float v2) throws StsException
    {
   		if(!growable)
        	throw new StsException("StsFloatArray2d.add(float,float,float)", "Can't use add method for fixed array.");

        if(nPntsAlloc == nPntsUsed) resize(nPntsUsed+incPnts);

        if(dimPnt < 3)
        {
        	throw new StsException("StsFloatArray2d.add(float,float,float)", "Point dimension: ",
            					   3, "> array dimension: ", dimPnt);
        }
        else
        {
        	array[nPntsUsed][0] = v0;
        	array[nPntsUsed][1] = v1;
        	array[nPntsUsed][2] = v2;
        	nPntsUsed++;
        }
    }

    public void setStsPointAt(int index, StsPoint point) throws StsException
    {
        setStsFloatsAt(index, point.v, 3);
    }

    public void setStsFloatsAt(int index, float[] v, int pointdimPnt) throws StsException
    {
        if(index >= nPntsAlloc)
        	throw new StsException("StsFloatArray2d.setStsPointAt(index,float[],int)", "Requested index: ", index, " >= number of points allocated: ", nPntsAlloc);

        if(dimPnt < pointdimPnt)
        	throw new StsException("StsFloatArray2d.setStsPoint3vAt(index,float[],int)", "Point dimension: ",
            					   pointdimPnt, "> array dimension: ", dimPnt);
        else
        {
        	for(int n = 0; n < pointdimPnt; n++)
        		array[index][n] = v[n];
        }
    }

    public void resize(int newNPnts)
    {
    	if(newNPnts <= nPntsAlloc) return;

        if(newNPnts < nPntsAlloc + incPnts)
        	newNPnts = nPntsAlloc + incPnts;

        float[][] oldArray = array;
        array = new float[newNPnts][dimPnt];

        for(int n = 0; n < nPntsAlloc; n++)
        	for(int m = 0; m < dimPnt; m++)
            	array[n][m] = oldArray[n][m];

//      System.out.println("StsFloatArray2d, resize from " + nPntsAlloc + " to " + newNPnts);

      nPntsAlloc = newNPnts;
    }

   	public int getNumberOfPointsUsed()
   	{
    	return nPntsUsed;
   	}

   	public void setNumberOfPointsUsed(int nPntsUsed)
   	{
    	this.nPntsUsed = nPntsUsed;
   	}

    public void setInterpolatedPoint(int index, float[] p0, float[] p1, float f) throws StsException
    {
    	int length0 = p0.length;
    	int length1 = p1.length;

        if(length0 != length1)
        	throw new StsException("StsFloatArray2d.setInterpolatedPoint()", "Points to be interpolated have different lengths: ", length0, " and ", length1);

        if(length0 <= 0)
       		throw new StsException("StsFloatArray2d.setInterpolatedPoint()", "Points to be interpolated have zero lengths: ", length0);

        if(length0 > dimPnt)
       	throw new StsException("StsFloatArray2d.setInterpolatedPoint()", "Point dimension: ", length0, "> array dimension: ", dimPnt);

    	for(int n = 0; n < length0; n++)
        	array[index][n] = p0[n] + f*(p1[n] - p0[n]);
    }

    public float normalizeReturnLength(int index)
   	{
    	float[] v = array[index];

    	double sumSq = (double)(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    	float vLength = (float)Math.sqrt(sumSq);

        if(vLength > 0.0)
        {
        	v[0] /= vLength;
        	v[1] /= vLength;
        	v[2] /= vLength;
        }

        return vLength;
    }

	public static void main(String[] args)
	{
    	StsFloatArray2d StsFloatArray;

    	try
        {
			StsFloatArray = new StsFloatArray2d(10, 10, 3);
        }
        catch (StsException e)
        {
        	System.out.println("StsException: " + e.getMessage());
            return;
        }

       	float x = 0.1f;
        float y = 0.2f;
        float z = 0.3f;

        for(int n = 0; n < 10; n++)
        {
        	try
            {
        		StsFloatArray.add(new float[] {x, y, z} );
            }
         	catch (StsException e)
            {
        		System.out.println("StsException: " + e.getMessage());
                return;
            }

            x += 1.0f; y += 1.0f; z += 1.0f;
        }

        for(int n = 10; n < 20; n++)
        {
        	try
            {
        		StsFloatArray.add(x, y, z);
            }
         	catch (StsException e)
            {
        		System.out.println("StsException: " + e.getMessage());
                return;
            }

            x += 1.0f; y += 1.0f; z += 1.0f;
        }

        try
        {
        	StsFloatArray.setStsPointAt(10, new StsPoint(110.1f, 110.2f, 110.3f) );
        }
        catch (StsException e)
        {
        	System.out.println("StsException: " + e.getMessage());
            return;
        }

      	float[][] array = StsFloatArray.array;

//	Note that this doesn't work:
//      StsPoint[] points = (StsPoint[])vector;

		int size = StsFloatArray.getNumberOfPointsUsed();
        System.out.println("number of points used: " + size);

        for(int n = 0; n < 20; n++)
        	System.out.println(n + " " + array[n][0] + " " + array[n][1] + " " + array[n][2]);
	}
}
