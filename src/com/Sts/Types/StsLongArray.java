
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

public class StsLongArray
{
	private boolean growable;		// this is growable array (not fixed)
    public int nPntsUsed;			// no of points currently used
	public int nPntsAlloc;			// number of allocated points
    public int dimPnt;              // dimension of points
    private int incPnts;            // 0 if fixed; number to increment by if growable
    public long array[] = null;  // allocated 2-dimensional array

	public StsLongArray(int nPnts, int incPnts)
	{
    	growable = (incPnts > 0);

        if(growable)
        	nPntsUsed = 0;
        else
        	nPntsUsed = nPnts;

    	this.incPnts = incPnts;

        try
        {
        	array = new long[nPnts];
            nPntsAlloc = nPnts;
        }
        catch (Exception e)
    	{
    		StsException.outputWarningException(this, "constructor", "Error allocating 2d float array of size : " + nPnts, e );
        }
	}

	public StsLongArray(int nPnts) throws StsException
	{
    	// This constructTraceAnalyzer makes a fixed array

    	growable = false;
       	nPntsUsed = nPnts;
        incPnts = 0;
        nPntsUsed = nPnts;
        this.dimPnt = dimPnt;

        try
        {
        	array = new long[nPnts];
        	nPntsAlloc = nPnts;
        }
        catch (Exception e)
    	{
    		throw new StsException("StsFloatArray2d.constructor1", "Can't allocate 2d float array of size : ", nPnts, " by ", dimPnt );
        }
	}

    /** accessors */
    public void setValues(long[] values) { array = values; }
    public long[] getValues() { return array; }

    public long elementAt(int index) throws StsException
    {
        if(index >= nPntsUsed)
        	throw new StsException("StsLongArray.elementAt", "Requested index: ", index, ">= no current points: ", nPntsUsed);
        else
        	return array[index];
    }

    public void add(long point) throws StsException
    {
    	if(!growable)
        	throw new StsException("StsFloatArray2d.add(float[])", "Can't use add method for fixed array.");
        if(nPntsAlloc == nPntsUsed) resize(nPntsAlloc+incPnts);
        array[nPntsUsed] = point;
        nPntsUsed++;
    }

    public void resize(int newNPnts)
    {
    	if(newNPnts <= nPntsAlloc) return;

        if(newNPnts < nPntsAlloc + incPnts)
        	newNPnts = nPntsAlloc + incPnts;

        long[] oldArray = array;
        array = new long[newNPnts];

        for(int n = 0; n < nPntsAlloc; n++)
            array[n] = oldArray[n];

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

    public long[] trimmedArray()
    {
        long[] newArray = new long[nPntsUsed];
        System.arraycopy(array, 0, newArray, 0, nPntsUsed);
        return newArray;
    }
}