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

public class StsFloatVector
{
    // constants
    static public final int MONOTONIC_UNKNOWN = 0;
    static public final int MONOTONIC_NOT = 1;
    static public final int MONOTONIC_INCR = 2;
    static public final int MONOTONIC_DECR = 3;
    static public final float largeFloat = Float.MAX_VALUE;
    // value range
    private float maxValue = -Float.MAX_VALUE;
    private float minValue = Float.MAX_VALUE;
    private float[] vector;
    private int size;
    private int growIncr = 1;
    private int capacity;
    private int monotonic = MONOTONIC_UNKNOWN;

	private float nullValue = StsParameters.nullValue;

    /** constructor for growable vector */
    public StsFloatVector()
    {
    }

    public StsFloatVector(int length, int incr)
    {
        vector = (length == 0) ? null : new float[length];
        size = 0;
        capacity = length;
        growIncr = incr;
    }

    /** constructor for growable vector with initial values to set */
    public StsFloatVector(float[] values, int length, int incr)
    {
        setValues(values, length);
        capacity = (values == null) ? length : values.length;
        growIncr = incr;
    }

    /** constructor for fixed length vector */
    public StsFloatVector(int length)
    {
        this(length, 0);
    }

    /** constructor for fixed length vector with initial values to set */
    public StsFloatVector(float[] values)
    {
        this(values, values.length, 0);
    }

    /** constructor for evenly-spaced vector from min to max by inc */
    public StsFloatVector(int nValues, float firstValue, float lastValue, float incValue)
    {
        vector = new float[nValues];
        float value = firstValue;
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

    /** set an array */
    public void setValues(float[] values, int length)
    {
        vector = values;
        size = (vector == null) ? 0 : Math.min(vector.length, length);
    }

    public void setValues(float[] values)
    {
        int length = (values == null) ? 0 : values.length;
        setValues(values, length);
    }

    public void resetVector()
    {
        size = 0;
        maxValue = -Float.MAX_VALUE;
        minValue = Float.MAX_VALUE;
        vector = new float[1];
        capacity = 1;
    }
    /* return the values array */

    public float[] getValues()
    {
        if(vector == null) return new float[0];
        else return vector; 
    }

    public void deleteValues() { vector = null; size = 0; }

    /* get value by index */

    public float getElement(int index)
    {
        try { return vector[index]; }
        catch (Exception e) { return StsParameters.nullValue; }
    }

    // set value by index with optional min/max & monotonic calculations
    // assume value is not null
    public boolean setElement(int index, float value)
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
    
    /** how much to grow array on append? */
    public void setGrowIncrement(int inc) { growIncr = inc; }
    
    /** append value to end of growable vector with optional min/max & monotonic calculation */
    public boolean append(float value)
    {
        if (size >= capacity)
        {
            if (growIncr == 0) return false;
            if(capacity == 0)
            	capacity = vector.length;            
            capacity += growIncr;
            float[] oldVector = vector;
            vector = new float[capacity];
            System.arraycopy(oldVector, 0, vector, 0, size);
        }
        vector[size++] = value;
        return true;
    }

    public float getMinValue() { return minValue; }
    public float getMaxValue() { return maxValue; }
    public void setMinValue(float minValue) { this.minValue = minValue; }
    public void setMaxValue(float maxValue) { this.maxValue = maxValue; }

    public void adjustRange(float adjustment)
    {
        if(minValue != largeFloat) minValue += adjustment;
        if(maxValue != -largeFloat) maxValue += adjustment;
    }

    public float getIncrement()
    {
        if (vector == null) return largeFloat;
        return vector[1] - vector[0];
    }

    /* calculate min and max for all the values */

    public boolean setMinMax()
    {
        if (vector == null || size == 0) return false;

        minValue = largeFloat;
        maxValue = -largeFloat;
        for (int i = 0; i < size; i++)
        {
            float value = vector[i];
            if(value != StsParameters.nullValue)
			{
                if (value < minValue) minValue = value;
                if (value > maxValue) maxValue = value;
            }
        }
        return true;
    }

	/** Sets the min-max range for the vector and
     *  converts values which equal the given nullValue to the standard nullValue
     *  and checks whether the complete vector is null.
     *  @param nullValue compare values to this nullValue
	 *  @return true if all values are null
	 */
	public boolean setMinMaxAndNulls(float nullValue)
	{
		if (vector == null || size == 0) return false;

		boolean allNull = true;
		minValue = largeFloat;
		maxValue = -largeFloat;
		for (int i = 0; i < size; i++)
		{
			float value = vector[i];
			if(value <= -largeFloat || value >= largeFloat)
				vector[i] = StsParameters.nullValue;
			else if(value == nullValue)
                vector[i] = StsParameters.nullValue;
            else
			{
				allNull = false;
				if (value < minValue) minValue = value;
				if (value > maxValue) maxValue = value;
			}
		}
		return allNull;
	}

    /** calculate monotonic value for values */
    public void checkMonotonic()
    {
        if (size <= 1)
        {
            monotonic = MONOTONIC_UNKNOWN;
            return;
        }
        if (vector[0] == vector[1]) // can't be monotonic
        {
            monotonic = MONOTONIC_NOT;
            return;
        }
        monotonic = (vector[1] > vector[0]) ? MONOTONIC_INCR : MONOTONIC_DECR;
        // see if increments stay in same direction
        float lastVal = vector[1];
        for (int i = 2; i < size; i++)
        {
            float val = vector[i];
            checkMonotonic(lastVal, val);
            lastVal = val;
        }
        // if we get to here, the values are monotonic
    }

    private void checkMonotonic(float lastVal, float val)
    {
        if (val < lastVal)
        {
            if (monotonic == MONOTONIC_INCR) monotonic = MONOTONIC_NOT;
            return;
        }
        else if (val > lastVal)
        {
            if (monotonic == MONOTONIC_DECR) monotonic = MONOTONIC_NOT;
            return;
        }
        else
            monotonic = MONOTONIC_NOT;
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
        if (monotonic == MONOTONIC_UNKNOWN) checkMonotonic();
        return monotonic != MONOTONIC_NOT;
    }

    /** get min/max indices enclosed by a range of values */
    public int[] getIndicesInValueRange(float minValue, float maxValue)
    {
        if (!isMonotonic()) return null; // range checks not valid
        if (size == 0) return null;
        float sign = (monotonic == MONOTONIC_DECR) ? -1.0f : 1.0f;
        if (sign < 0) { minValue *= sign; maxValue *= sign; }
        if (minValue > maxValue)
        {
            float temp = minValue;
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
    
    public void trimToSize()
    {
        float[] newVector = new float[size];
        System.arraycopy(vector, 0, newVector, 0, size);
        vector = newVector;
        capacity = size;
    }

    /** set/get null value */
    //    public void setNullValue(float nullValue) { this.nullValue = nullValue; }
    //    public float getNullValue() { return nullValue; }
    public static void main(String[] args)
    {
        float[] vector = { 2.0f, 3.0f, 4.0f, 4.1f, 4.3f, 5.1f, 11.0f };
        StsFloatVector fv = new StsFloatVector(vector);
        float minValue = fv.getMinValue();
        float maxValue = fv.getMaxValue();
        int monotonic = fv.getMonotonic();
        boolean isMonotonic = fv.isMonotonic();
        float[] range = { 4.0f, 5.3f };
        int[] indexRange = fv.getIndicesInValueRange(range[0], range[1]);
        System.out.println("\nStsFloatVector test:\n");
        System.out.print("\tvector: ");
        for (int i = 0; i < vector.length; i++) System.out.print(vector[i] + " ");
        System.out.print(" \n");
        System.out.println("\tminValue = " + minValue);
        System.out.println("\tmaxValue = " + maxValue);
        System.out.println("\tmonotonic = " + monotonic);
        System.out.println("\tisMonotonic = " + isMonotonic);
        System.out.print("\trange: " + range[0] + " - " + range[1]);
        if (indexRange != null) System.out.println(",  index range: " + indexRange[0] + " - " + indexRange[1]);
        else
            System.out.println(", index range: null");
        range[0] = 1.0f; range[1] = 4.3f;
        indexRange = fv.getIndicesInValueRange(range[0], range[1]);
        System.out.print("\trange: " + range[0] + " - " + range[1]);
        if (indexRange != null) System.out.println(",  index range: " + indexRange[0] + " - " + indexRange[1]);
        else
            System.out.println(", index range: null");
        range[0] = 13.0f; range[1] = 4.3f;
        indexRange = fv.getIndicesInValueRange(range[0], range[1]);
        System.out.print("\trange: " + range[0] + " - " + range[1]);
        if (indexRange != null) System.out.println(",  index range: " + indexRange[0] + " - " + indexRange[1]);
        else
            System.out.println(", index range: null");
        range[0] = 13.0f; range[1] = 14.3f;
        indexRange = fv.getIndicesInValueRange(range[0], range[1]);
        System.out.print("\trange: " + range[0] + " - " + range[1]);
        if (indexRange != null) System.out.println(",  index range: " + indexRange[0] + " - " + indexRange[1]);
        else
            System.out.println(", index range: null");
        range[0] = 10.0f; range[1] = 14.3f;
        indexRange = fv.getIndicesInValueRange(range[0], range[1]);
        System.out.print("\trange: " + range[0] + " - " + range[1]);
        if (indexRange != null) System.out.println(",  index range: " + indexRange[0] + " - " + indexRange[1]);
        else
            System.out.println(", index range: null");
        range[0] = 6.0f; range[1] = 7.0f;
        indexRange = fv.getIndicesInValueRange(range[0], range[1]);
        System.out.print("\trange: " + range[0] + " - " + range[1]);
        if (indexRange != null) System.out.println(",  index range: " + indexRange[0] + " - " + indexRange[1]);
        else
            System.out.println(", index range: null");
        float[] vector2 = { 11.0f, 5.1f, 4.3f, 4.1f, 4.0f, 3.0f, 2.0f };
        fv.setValues(vector2, vector2.length);
        minValue = fv.getMinValue();
        maxValue = fv.getMaxValue();
        monotonic = fv.getMonotonic();
        isMonotonic = fv.isMonotonic();
        range[0] = 4.0f; range[1] = 5.3f;
        indexRange = fv.getIndicesInValueRange(range[0], range[1]);
        System.out.println("\nStsFloatVector test2:\n");
        System.out.print("\tvector2: ");
        for (int i = 0; i < vector2.length; i++) System.out.print(vector2[i] + " ");
        System.out.print(" \n");
        System.out.println("\tminValue = " + minValue);
        System.out.println("\tmaxValue = " + maxValue);
        System.out.println("\tmonotonic = " + monotonic);
        System.out.println("\tisMonotonic = " + isMonotonic);
        System.out.print("\trange: " + range[0] + " - " + range[1]);
        if (indexRange != null) System.out.println(",  index range: " + indexRange[0] + " - " + indexRange[1]);
        else
            System.out.println(", index range: null");
        range[0] = 1.0f; range[1] = 4.3f;
        indexRange = fv.getIndicesInValueRange(range[0], range[1]);
        System.out.print("\trange: " + range[0] + " - " + range[1]);
        if (indexRange != null) System.out.println(",  index range: " + indexRange[0] + " - " + indexRange[1]);
        else
            System.out.println(", index range: null");
        range[0] = 13.0f; range[1] = 4.3f;
        indexRange = fv.getIndicesInValueRange(range[0], range[1]);
        System.out.print("\trange: " + range[0] + " - " + range[1]);
        if (indexRange != null) System.out.println(",  index range: " + indexRange[0] + " - " + indexRange[1]);
        else
            System.out.println(", index range: null");
        range[0] = 13.0f; range[1] = 14.3f;
        indexRange = fv.getIndicesInValueRange(range[0], range[1]);
        System.out.print("\trange: " + range[0] + " - " + range[1]);
        if (indexRange != null) System.out.println(",  index range: " + indexRange[0] + " - " + indexRange[1]);
        else
            System.out.println(", index range: null");
        range[0] = 10.0f; range[1] = 14.3f;
        indexRange = fv.getIndicesInValueRange(range[0], range[1]);
        System.out.print("\trange: " + range[0] + " - " + range[1]);
        if (indexRange != null) System.out.println(",  index range: " + indexRange[0] + " - " + indexRange[1]);
        else
            System.out.println(", index range: null");
        range[0] = 6.0f; range[1] = 7.0f;
        indexRange = fv.getIndicesInValueRange(range[0], range[1]);
        System.out.print("\trange: " + range[0] + " - " + range[1]);
        if (indexRange != null) System.out.println(",  index range: " + indexRange[0] + " - " + indexRange[1]);
        else
            System.out.println(", index range: null");
    }
}
