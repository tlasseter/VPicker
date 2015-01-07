package com.Sts.IO;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

import com.Sts.Utilities.*;

import java.io.*;

public class StsMappedByteBuffer extends StsMappedBuffer
{
//    public ByteBuffer buffer = null;
	public byte[] scratchBytes = null;
	float dataMin, dataMax;
    float scale, scaleOffset;

    private StsMappedByteBuffer(String directory, String filename, String mode, float dataMin, float dataMax) throws FileNotFoundException
    {
		super(directory, filename, mode);
        this.dataMin = dataMin;
		this.dataMax = dataMax;
        scale = StsMath.unsignedIntToFloatScale(dataMin, dataMax);
        scaleOffset = dataMin;
    }

    private StsMappedByteBuffer(String directory, String filename, String mode) throws FileNotFoundException
    {
		super(directory, filename, mode);
        scale = StsMath.unsignedIntToFloatScale(dataMin, dataMax);
        scaleOffset = dataMin;
    }

    static public StsMappedByteBuffer constructor(String directory, String filename, String mode, float dataMin, float dataMax)
	{
		try
		{
			return new StsMappedByteBuffer(directory, filename, mode, dataMin, dataMax);
        }
		catch(Exception e)
		{
			StsException.outputException("StsMappedByteBuffer.constructor() failed.", e, StsException.WARNING);
			return null;
		}
    }

	static public StsMappedByteBuffer constructor(String directory, String filename, String mode)
	{
		try
		{
			return new StsMappedByteBuffer(directory, filename, mode);
        }
		catch(Exception e)
		{
			StsException.outputException("StsMappedByteBuffer.constructor() failed.", e, StsException.WARNING);
			return null;
		}
    }

    public boolean map(long position, long nSamples)
	{
		try
		{
            // for consistence with float and double mappedByteBuffers, we first allocate the base type: byte buffer;
            // for bytes, floats, and doubles, we then create a buffer of that specific type; here it is redundant since
            // the base buffer is bytes and the desired buffer is bytes.
            byteBuffer = channel.map(mapMode, position, nSamples);
//            buffer = byteBuffer.asReadOnlyBuffer();
            return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsMappedFloatBuffer.map() failed for position: " + position + " nSamples: " + nSamples,
										 e, StsException.WARNING);
			return false;
		}
	}

    public long getBufferPosition()
    {
        if(byteBuffer == null) return 0;
        return byteBuffer.position();
    }

	final public long getCapacity()
	{
        if(byteBuffer == null) return 0;
        return byteBuffer.capacity();
	}

	final public void put(byte[] bytes)
	{
		byteBuffer.put(bytes);
	}

	final public void put(byte b)
	{
		byteBuffer.put(b);
	}

	final public void get(byte[] bytes)
	{
		byteBuffer.get(bytes);
	}

	final public void get(float[] floats)
	{
		int nValues = floats.length;
		if(scratchBytes == null || scratchBytes.length != nValues)
			scratchBytes = new byte[nValues];
		byteBuffer.get(scratchBytes);
		for(int n = 0; n < nValues; n++)
		{
			if (scratchBytes[n] == -1)
				floats[n] = StsParameters.nullValue;
			else
            {
                int unsignedInt = StsMath.signedByteToUnsignedInt(scratchBytes[n]);
                floats[n] = StsMath.unsignedIntToFloat(unsignedInt, scale, scaleOffset);
            }
        }
    }
      /** if the mappedBuffer capacity is exceeded or doesn't exist, it is cleared and remapped so writing can continue */
    public void checkPut(byte[] bytes, long blockSize)
	{
		try
		{
            int nBytes = bytes.length;
            checkPosition(nBytes, blockSize);
            byteBuffer.put(bytes);
		}
		catch(Exception e)
		{
			StsException.systemError(this, "put(float[])", "Position: " + byteBuffer.position() + " buffer remaining: " + byteBuffer.remaining() + " array length " + bytes.length);
		}
	}

    /** if the mappedBuffer capacity is exceeded or doesn't exist, it is cleared and remapped so reading can continue */
    final public void checkGet(byte[] bytes, long blockSize)
	{
        int nBytes = bytes.length;
        checkPosition(nBytes, blockSize);
        byteBuffer.get(bytes);
	}

	final public float getFloat()
	{
		byte b = byteBuffer.get();
		if(b == -1)
            return StsParameters.nullValue;
		else
            return StsMath.signedByteToFloatWithScale(b, scale, scaleOffset);
	}

	final public void get(double[] doubles)
	{
		int nValues = doubles.length;
		if(scratchBytes == null || scratchBytes.length != nValues)
			scratchBytes = new byte[nValues];
		byteBuffer.get(scratchBytes);
		for(int n = 0; n < nValues; n++)
		{
			if (scratchBytes[n] == -1)
            {
                int prevGoodIndex = n-1;
                int nextGoodIndex = nValues;
                n++;
                for(; n < nValues; n++)
                {
                    if(scratchBytes[n] != -1)
                    {
                        nextGoodIndex = n;
                        break;
                    }
                }
                if(prevGoodIndex < 0)
                {
                    if(nextGoodIndex == nValues) // no non-nulls in trace; set values to 0
                    {
                        for(int nn = 0; nn < nValues; nn++)
                            doubles[nn] = 0.0;
                    }
                    else // trace begins with nulls; set values to first non-null value
                    {
                        double value = getScaledDoubleValue(scratchBytes[nextGoodIndex]);
                        for(int nn = 0; nn <= nextGoodIndex; nn++)
                            doubles[nn] = value;
                    }
                }
                else
                {
                    double prevValue = getScaledDoubleValue(scratchBytes[prevGoodIndex]);
                    if(nextGoodIndex == nValues) // trace ends with nulls; set values to 0
                    {
                        for(int nn = prevGoodIndex+1; nn < nValues; nn++)
                            doubles[nn] = prevValue;
                    }
                    else // trace has nulls in this section; interpolate values
                    {
                        double nextValue = getScaledDoubleValue(scratchBytes[nextGoodIndex]);
                        double dValue = (nextValue - prevValue)/(nextGoodIndex - prevGoodIndex);
                        double value = prevValue + dValue;
                        for(int nn = prevGoodIndex+1; nn < nextGoodIndex; nn++, value += dValue)
                            doubles[nn] = value;
                        doubles[nextGoodIndex] = nextValue;
                    }
                }
            }
			else
                doubles[n] = getScaledDoubleValue(scratchBytes[n]);
		}
	}

	final public double getScaledDoubleValue(byte b)
	{
		double scaledValue = (double)StsMath.signedByteToUnsignedInt(b);
		return dataMin + (scaledValue / 254) * (dataMax - dataMin);
	}

	final public double getDouble()
	{
		byte b = byteBuffer.get();
		if(b == -1) return StsParameters.nullValue;
		else       return getScaledDoubleValue(b);
	}

	final public void position(int position)
	{
		byteBuffer.position(position);
	}

	final public void rewind() { byteBuffer.rewind(); }

	public void clear0()
	{
		if(byteBuffer != null) byteBuffer.clear();
	}

	public void force()
	{
		if(byteBuffer != null) byteBuffer.force();
	}
}
