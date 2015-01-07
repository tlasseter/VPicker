package com.Sts.IO;

import com.Sts.Utilities.*;

import java.io.*;
import java.nio.*;


public class StsMappedDoubleBuffer extends StsMappedBuffer
{
	public DoubleBuffer doubleBuffer = null;
	double[] scratchDoubles = null;

	private StsMappedDoubleBuffer(String directory, String filename, String mode) throws FileNotFoundException
	{
		super(directory, filename, mode);
	}

	static public StsMappedDoubleBuffer constructor(String directory, String filename, String mode) throws FileNotFoundException
	{
		return new StsMappedDoubleBuffer(directory, filename, mode);
	}

	public boolean map(long position, long nSamples)
	{
		try
		{
			this.position = position;
			byteBuffer = channel.map(mapMode, 8*position, 8*nSamples);
			doubleBuffer = byteBuffer.asDoubleBuffer();
			return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsMappedDoubleBuffer.map() failed for position: " + position + " nSamples: " + nSamples,
										 e, StsException.WARNING);
			return false;
		}
	}

    public long getBufferPosition()
    {
        if(doubleBuffer == null) return 0;
        return doubleBuffer.position();
    }

    final public long getCapacity()
	{
        if(doubleBuffer == null) return 0;
        return doubleBuffer.capacity();
	}
	public void put(double[] doubles, int offset, int length)
	{
		doubleBuffer.put(doubles, offset, length);
	}
	public void put(double d)
	{
		doubleBuffer.put(d);
	}

	final public void get(double[] doubles)
	{
		doubleBuffer.get(doubles);
	}

	final public double getDouble()
	{
		return doubleBuffer.get();
	}
	final public void get(float[] floats)
	{
		int nValues = floats.length;
		if(scratchDoubles == null || scratchDoubles.length != nValues)
			scratchDoubles = new double[nValues];
		doubleBuffer.get(scratchDoubles);
		for(int n = 0; n < nValues; n++)
		{
			if (scratchDoubles[n] == StsParameters.doubleNullValue) floats[n] = StsParameters.nullValue;
			else floats[n] = (float)scratchDoubles[n];
		}
	}

    /** if the mappedBuffer capacity is exceeded or doesn't exist, it is cleared and remapped so writing can continue */
    public void checkPut(double[] doubles, long blockSize)
	{
		try
		{
            int nDoubles = doubles.length;
            checkPosition(nDoubles, blockSize);
            doubleBuffer.put(doubles);
		}
		catch(Exception e)
		{
			StsException.systemError(this, "put(double[])", "Position: " + doubleBuffer.position() + " buffer remaining: " + doubleBuffer.remaining() + " array length " + doubles.length);
		}
	}

    /** if the mappedBuffer capacity is exceeded or doesn't exist, it is cleared and remapped so reading can continue */
    final public void checkGet(double[] doubles, long blockSize)
	{
        int nDoubles = doubles.length;
        checkPosition(nDoubles, blockSize);
        doubleBuffer.get(doubles);
	}
	final public float getFloat()
	{
		return (float) doubleBuffer.get();
	}

	final public void position(int position)
	{
		doubleBuffer.position(position);
	}

	final public void rewind() { doubleBuffer.rewind(); }

	public void clear0()
	{
		if(doubleBuffer != null) doubleBuffer.clear();
	}
}
