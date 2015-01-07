package com.Sts.Utilities;

/**
 * Title:        Sts Well Viewer Description:  Well Model-Viewer Copyright:    Copyright (c) 2001
 * Company:      4D Systems LLC
 * @author T.J.Lasseter
 * @version 1.0
 */
// Allocates a vector of float values
// These arrays can be either fixed or growable.
// If fixed, call them with incr of 0.  length is set
// and array is flagged as fixed. Values are set using set.
// If growable, call them with positive incr.  length is set
// and values are added to the array with the add method.

public class StsLongVector
{
    // constants
    static public final int MONOTONIC_UNKNOWN = 0;
    static public final int MONOTONIC_NOT = 1;
    static public final int MONOTONIC_INCR = 2;
    static public final int MONOTONIC_DECR = 3;
    static public final long largeLong = StsParameters.largeLong;
    // value range
    private long maxValue = -largeLong;
    private long minValue = largeLong;
    private long[] vector;
    private int size;
    private int growIncr = 1;
    private int capacity;
    private int monotonic = MONOTONIC_UNKNOWN;

	private float nullValue = StsParameters.nullValue;

    /** constructor for growable vector */
    public StsLongVector()
    {
    }

    public StsLongVector(int length, int incr)
    {
        vector = (length == 0) ? null : new long[length];
        size = 0;
        capacity = length;
        growIncr = incr;
    }

    /** constructor for growable vector with initial values to set */
    public StsLongVector(long[] values, int length, int incr)
    {
        setValues(values, length);
        capacity = (values == null) ? length : values.length;
        growIncr = incr;
    }

    /** constructor for fixed length vector */
    public StsLongVector(int length)
    {
        this(length, 0);
    }

    /** constructor for fixed length vector with initial values to set */
    public StsLongVector(long[] values)
    {
        this(values, values.length, 0);
    }

    /** constructor for evenly-spaced vector from min to max by inc */
    public StsLongVector(int nValues, long firstValue, long lastValue, long incValue)
    {
        vector = new long[nValues];
        long value = firstValue;
        for (int n = 0; n < nValues; n++, value += incValue)
            vector[n] = value;
        if (firstValue <= lastValue)
        {
            minValue = firstValue;
            maxValue = lastValue;
        }
        else
        {
            minValue = lastValue;
            maxValue = firstValue;
        }
    }

    /** how many values are in the vector? */
    public int getSize() { return size; }

    /** how much to grow array on append? */
    public void setGrowIncrement(int inc) { growIncr = inc; }

    /** set an array */
    public void setValues(long[] values, int length)
    {
        vector = values;
        size = (vector == null) ? 0 : Math.min(vector.length, length);
    }

    public void setValues(long[] values)
    {
        int length = (values == null) ? 0 : values.length;
        setValues(values, length);
    }

    /* return the values array */

    public long[] getValues() { return vector; }

    public float[] getValuesAsRelativeFloats(long startPoint)
    {
    	float[] flts = new float[vector.length];
    	for(int i=0; i<flts.length; i++)
    		flts[i] = (float)((vector[i] - startPoint)/1000);
    	return flts;
    }
    public void deleteValues() { vector = null; size = 0; }

    /* get value by index */

    public long getElement(int index)
    {
        try { return vector[index]; }
        catch (Exception e) { return StsParameters.nullLongValue; }
    }

    public long getFirst() { return vector[0]; }
    public long getLast() { return vector[size-1]; }

    // set value by index with optional min/max & monotonic calculations
    // assume value is not null
    public boolean setElement(int index, long value)
    {
        try
        {
            vector[index] = value;
            if (value < minValue) minValue = value;
            if (value > maxValue) maxValue = value;
            checkMonotonic(index);
            return true;
        }
        catch (Exception e) { return false; }
    }

    /** append value to end of growable vector with optional min/max & monotonic calculation */
    public boolean append(long value)
    {
        if (size >= capacity)
        {
            if (growIncr == 0) return false;
            if(capacity == 0)
            	capacity = vector.length;
            capacity += growIncr;
            long[] oldVector = vector;
            vector = new long[capacity];
            System.arraycopy(oldVector, 0, vector, 0, size);
        }
        vector[size++] = value;
        return true;
    }

    public long getMinValue() { return minValue; }
    public long getMaxValue() { return maxValue; }
    public void setMinValue(long minValue) { this.minValue = minValue; }
    public void setMaxValue(long maxValue) { this.maxValue = maxValue; }

    public void adjustRange(float adjustment)
    {
        if(minValue != largeLong) minValue += adjustment;
        if(maxValue != -largeLong) maxValue += adjustment;
    }

    public float getIncrement()
    {
        if (vector == null) return largeLong;
        return vector[1] - vector[0];
    }

    /* calculate min and max for all the values */

    public boolean setMinMax()
    {
        if (vector == null || size == 0) return false;

        minValue = largeLong;
        maxValue = -largeLong;
        for (int i = 0; i < size; i++)
        {
            long value = vector[i];
            if(value != StsParameters.nullValue)
			{
                if (value < minValue) minValue = value;
                if (value > maxValue) maxValue = value;
            }
        }
        return true;
    }

	/** @param nullValue compare values to this nullValue
	 *  @return true if all values are null
	 */
	public boolean setMinMaxAndNulls(long nullValue)
	{
		if (vector == null || size == 0) return false;

		boolean isNull = true;
		minValue = largeLong;
		maxValue = -largeLong;
		for (int i = 0; i < size; i++)
		{
			long value = vector[i];
			if(value <= -largeLong || value >= largeLong)
				vector[i] = nullValue;
			else if(value != nullValue)
			{
				isNull = false;
				if (value < minValue) minValue = value;
				if (value > maxValue) maxValue = value;
			}
		}
		return isNull;
	}

    public void resetVector()
    {
        size = 0;
        maxValue = -largeLong;
        minValue = largeLong;
        vector = new long[1];
        capacity = 1;
    }

    /** calculate if values are in monotonically increasing or decreasing order.
     *  Two successive values which are equal will be ignored.
     */
    public int checkMonotonic()
    {
        if(vector == null || size <= 1)
            monotonic = MONOTONIC_UNKNOWN;
        else
        {
            // see if increments stay in same direction
            long lastVal = vector[0];
            for (int i = 1; i < size; i++)
            {
                long val = vector[i];
                if(val != lastVal)
                {
                    checkMonotonic(lastVal, val);
                    if(monotonic == MONOTONIC_NOT)
                        return monotonic;
                    lastVal = val;
                }
            }
        }
        return monotonic;
    }

    private void checkMonotonic(long lastVal, long val)
    {
        if (val < lastVal)
        {
            if (monotonic == MONOTONIC_INCR)
                monotonic = MONOTONIC_NOT;
            else if(monotonic == MONOTONIC_UNKNOWN)
                monotonic = MONOTONIC_DECR;
        }
        else if (val > lastVal)
        {
            if (monotonic == MONOTONIC_DECR)
                monotonic = MONOTONIC_NOT;
            else if(monotonic == MONOTONIC_UNKNOWN)
                monotonic = MONOTONIC_INCR;
        }
    }

    private void checkMonotonic(int index)
    {
        if (index > 0) checkMonotonic(vector[index - 1], vector[index]);
        if (index < size - 1) checkMonotonic(vector[index], vector[index + 1]);
    }

    /** get monotonic value */
    public int getMonotonic() { return monotonic; }

    /** get monotonic status */
    public boolean isMonotonic()
    {
        if(monotonic == MONOTONIC_UNKNOWN) checkMonotonic();
        if(monotonic == MONOTONIC_NOT)
            System.out.println("Monotinically NOT");
        return monotonic != MONOTONIC_NOT;
    }

    /** get min/max indices enclosed by a range of values */
    public int[] getIndicesInValueRange(long minValue, long maxValue)
    {
        //if (!isMonotonic()) return null; // range checks not valid
        if (size == 0) return null;
        long sign = (monotonic == MONOTONIC_DECR) ? -1l : 1l;
        if (sign < 0) { minValue *= sign; maxValue *= sign; }
        if (minValue > maxValue)
        {
            long temp = minValue;
            minValue = maxValue;
            maxValue = temp;
        }
        int[] indexRange = new int[2];
        // get start index
        if (minValue < vector[0] * sign) indexRange[0] = 0;
        else if (minValue > vector[size - 1] * sign) indexRange[0] = size;
        else
        {
            for (int i = 0; i < size; i++)
            {
                if (minValue <= vector[i] * sign)
                {
                    indexRange[0] = i;
                    break;
                }
            }
        }
        // get end index
        if (maxValue < vector[0] * sign) indexRange[1] = -1;
        else if (maxValue > vector[size - 1] * sign) indexRange[1] = size - 1;
        else
        {
            for (int i = 0; i < size; i++)
            {
                if (maxValue == vector[i] * sign)
                {
                    indexRange[1] = i;
                    break;
                }
                else if (maxValue < vector[i] * sign)
                {
                    indexRange[1] = i - 1;
                    break;
                }
            }
        }
        if (indexRange[1] < indexRange[0]) return null;
        return indexRange;
    }

    public long[] trimToSize()
    {
        long[] newVector = new long[size];
        System.arraycopy(vector, 0, newVector, 0, size);
        vector = newVector;
        capacity = size;
        return vector;
    }

}
