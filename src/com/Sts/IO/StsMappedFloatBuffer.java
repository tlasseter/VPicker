package com.Sts.IO;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.nio.*;

public class StsMappedFloatBuffer extends StsMappedBuffer
{
    public FloatBuffer floatBuffer = null;
	public float[] scratchFloats = null;
    float dataMin, dataMax;

	private StsMappedFloatBuffer(String directory, String filename, String mode) throws FileNotFoundException
	{
		super(directory, filename, mode);
    }

    static public StsMappedFloatBuffer openRead(String directory, String filename)
	{
		try
		{
			return new StsMappedFloatBuffer(directory, filename, READ);
		}
		catch(Exception e)
		{
            StsException.outputWarningException(StsMappedFloatBuffer.class, "openRead", e);
			return null;
		}
	}

	static public StsMappedFloatBuffer openReadWrite(String directory, String filename)
	{
		try
		{
			return new StsMappedFloatBuffer(directory, filename, READ_WRITE);
		}
		catch(Exception e)
		{
            StsException.outputWarningException(StsMappedFloatBuffer.class, "openReadWrite", e);
			return null;
		}
	}

	static public StsMappedFloatBuffer constructor(String directory, String filename, String mode)
	{
		try
		{
			return new StsMappedFloatBuffer(directory, filename, mode);
		}
		catch(Exception e)
		{
//			StsException.outputException("StsMappedFloatBuffer.constructor() failed.", e, StsException.WARNING);
			return null;
		}
	}


    public boolean map(long position, long nSamples)
	{
        try
        {
            this.position = position;
            byteBuffer = channel.map(mapMode, 4*position, 4*nSamples);
            floatBuffer = byteBuffer.asFloatBuffer();
            return true;
        }
        catch(IOException ioe)
        {
            return false;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "map(position, nSamples)", e);
            return false;
        }
    }

    public boolean map(long position, long nSamples, StsProgressPanel progressPanel)
	{
        try
        {
            this.position = position;
            byteBuffer = channel.map(mapMode, 4*position, 4*nSamples);
            floatBuffer = byteBuffer.asFloatBuffer();
            return true;
        }
        catch(IOException ioe)
        {
            progressPanel.appendLine("    IO Exception: " + ioe.getMessage());
            return false;
        }
        catch(Exception e)
        {
            progressPanel.appendLine("    Exception: " + e.getMessage());
            return false;
        }
    }
    
    public long getBufferPosition()
    {
        if(floatBuffer == null) return 0;
        return floatBuffer.position();
    }


    final public long getCapacity()
	{
        if(floatBuffer == null) return 0;
        return floatBuffer.capacity();
	}

    public void put(float[] floats)
	{
		try
		{
			floatBuffer.put(floats);
		}
		catch(Exception e)
		{
			StsException.systemError(this, "put(float[])", "Position: " + floatBuffer.position() + " buffer remaining: " + floatBuffer.remaining() + " array length " + floats.length);
		}
	}

	public void put(float[] floats, int offset, int length)
	{
		floatBuffer.put(floats, offset, length);
	}

	public void put(float f)
	{
		floatBuffer.put(f);
	}

    final public void get(float[] floats)
	{
		floatBuffer.get(floats);
	}

    /** if the mappedBuffer capacity is exceeded or doesn't exist, it is cleared and remapped so writing can continue */
    public void checkPut(float[] floats, long blockSize)
	{
		try
		{
            int nFloats = floats.length;
            checkPosition(nFloats, blockSize);
            floatBuffer.put(floats);
            position += nFloats;
        }
		catch(Exception e)
		{
			StsException.systemError(this, "put(float[])", "Position: " + floatBuffer.position() + " buffer remaining: " + floatBuffer.remaining() + " array length " + floats.length);
		}
	}

    /** if the mappedBuffer capacity is exceeded or doesn't exist, it is cleared and remapped so writing can continue */
    public void checkPutDebug(float[] floats, long blockSize, String debugString)
	{
		try
		{
            int nFloats = floats.length;
            checkPositionDebug(nFloats, blockSize, debugString);
            floatBuffer.put(floats);
            position += nFloats;
        }
		catch(Exception e)
		{
			StsException.systemError(this, "put(float[])", "Position: " + floatBuffer.position() + " buffer remaining: " + floatBuffer.remaining() + " array length " + floats.length);
		}
	}

    /** if the mappedBuffer capacity is exceeded or doesn't exist, it is cleared and remapped so reading can continue */
    final public void checkGet(float[] floats, long blockSize)
	{
        int nFloats = floats.length;
        checkPosition(nFloats, blockSize);
        floatBuffer.get(floats);
        position += nFloats;
    }

    final public void checkGetDebug(float[] floats, long blockSize, String debugString)
	{
        int nFloats = floats.length;
        checkPositionDebug(nFloats, blockSize, debugString);
        floatBuffer.get(floats);
        position += nFloats;
    }
	final public float getFloat()
	{
		return floatBuffer.get();
	}

	final public void get(double[] doubles)
	{
		int nValues = doubles.length;
		if(scratchFloats == null || scratchFloats.length != nValues)
			scratchFloats = new float[nValues];
		floatBuffer.get(scratchFloats);
		for(int n = 0; n < nValues; n++)
		{
			if (scratchFloats[n] == StsParameters.nullValue) doubles[n] = StsParameters.nullValue;
			else doubles[n] = (double)scratchFloats[n];
		}
	}

	final public ByteBuffer getScaledByteBuffer(float dataMin, float dataMax, int nValues)
	{
		if(scratchFloats == null || scratchFloats.length != nValues)
			scratchFloats = new float[nValues];
		floatBuffer.get(scratchFloats);
        byte[] bytes =  StsMath.floatsToUnsignedBytes254(scratchFloats, dataMin, dataMax);
        return ByteBuffer.wrap(bytes);
    }

    final public double getDouble()
	{
		return (double) floatBuffer.get();
	}

	final public void position(int position)
	{
		floatBuffer.position(position);
	}

	final public void rewind() { floatBuffer.rewind(); }

	public void clear0()
	{
		if(floatBuffer != null) floatBuffer.clear();
	}
    static public void main(String[] args)
    {
        int nFloats = 101;
        int nDataBlocks = 10;
        long blockSize = 200;
        float[] floats = new float[nFloats];

        String directory = "c:\\";
        String filename = "test";
        StsMappedFloatBuffer floatBuffer = StsMappedFloatBuffer.openReadWrite(directory, filename);
//        floatBuffer.map(0, blockSize);
        int nn = 0;
        for(int i = 0; i < nDataBlocks; i++)
        {
            for(int n = 0; n < nFloats; n++, nn++)
                floats[n] = nn;
            floatBuffer.checkPutDebug(floats, blockSize, "output block " + i + ".");
        }
        floatBuffer.close();

        floatBuffer = StsMappedFloatBuffer.openRead(directory, filename);
//        floatBuffer.map(0, blockSize);
        for(int i = 0; i < nDataBlocks; i++)
        {
            floatBuffer.checkGetDebug(floats, blockSize, "input block " + i + ".");
            for(int n = 0; n < nFloats; n++)
                System.out.print(floats[n] + " ");
            System.out.println();
        }
    }
}
