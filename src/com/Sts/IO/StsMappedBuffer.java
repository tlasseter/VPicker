package com.Sts.IO;

import com.Sts.MVC.Main;
import com.Sts.Utilities.*;

import java.io.*;
import java.nio.*;
import java.nio.channels.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

abstract public class StsMappedBuffer
{
	public RandomAccessFile randomAccessFile = null;
	public MappedByteBuffer byteBuffer = null;
	public FileChannel channel;
	public FileChannel.MapMode mapMode = FileChannel.MapMode.READ_ONLY;
    public String pathname;
    public long position; // position in file in bytes

    static public final String READ = "r";
    static public final String WRITE = "w";
    static public final String READ_WRITE = "rw";
	public StsMappedBuffer(String directory, String filename, String mode) throws FileNotFoundException, IllegalArgumentException
	{
		try
		{
            if(Main.isJarDB)
            {
                String[] filenames = Main.jar.getFilenames();
                String rootDir = System.getProperty("user.home") + File.separator + "S2SCache" + File.separator;
                for(int i=0; i<filenames.length; i++)
                {
                    if(filenames[i].endsWith(filename))
                    {
                        File file = Main.jar.uncompressTo(rootDir, filename);
                        pathname = file.getPath();
                        break;
                    }
                }
            }
            else
            	pathname = directory + filename;
			randomAccessFile = new RandomAccessFile(pathname, mode);
			if (!mode.equals(READ)) mapMode = FileChannel.MapMode.READ_WRITE;
			channel = randomAccessFile.getChannel();
		}
		catch(Exception e)
		{
			StsException.systemError("StsMappedBuffer.constructor() failed. Illegal mode: " + mode + " must be \"r\" , \"w\" or \"rw\" ");
		}
	}

	public void clean()
	{
		StsToolkit.clean(byteBuffer);
	}

    public void clear()
	{
		try
		{
            if(byteBuffer == null) return;
            byteBuffer.force();
            StsToolkit.clean(byteBuffer);
            channel.force(false);
		}
		catch(Exception e)
		{
		}
	}
    public void clearDebug(String message, StsTimer timer)
	{
		try
		{
            timer.start();
            byteBuffer.force();
            timer.stopPrint(message + " byteBuffer.force()");

            timer.start();
            StsToolkit.clean(byteBuffer);
            timer.stopPrint(message + " clean(byteBuffer)");

            timer.start();
            channel.force(false);
            timer.stopPrint(message + " channel.force()");
        }
		catch(Exception e)
		{
		}
	}

    public void close()
	{
		try
		{
            clear();
			channel.close();
			randomAccessFile.close();
			randomAccessFile = null;
//			clear0();
		}
		catch(Exception e)
		{
		}
	}

    public void force()
    {
        if(byteBuffer != null) byteBuffer.force();
    }

    abstract public void rewind();
	abstract public void clear0();
	abstract public boolean map(long position, long nSamples);
	abstract public void get(float[] floats);
	abstract public float getFloat();
	abstract public void position(int position);
    abstract public long getBufferPosition();
	abstract public void get(double[] doubles);
	abstract public double getDouble();
    abstract public long getCapacity();

    public long getFilePosition() { return position; }

    protected boolean checkPosition(int nSamples, long blockSize)
    {
        long currentBufferPosition = getBufferPosition();
        long capacity = getCapacity();
        if(currentBufferPosition + nSamples > capacity)
        {
            long currentFilePosition = getFilePosition();
            clear();
            map(currentFilePosition, blockSize);
        }
        return true;
    }

    protected boolean checkPositionDebug(int nSamples, long blockSize, String debugString)
    {
        long currentBufferPosition = getBufferPosition();
        long currentFilePosition = getFilePosition();
        long capacity = getCapacity();
        if(currentBufferPosition + nSamples > capacity)
        {
            clear();
            map(currentFilePosition, blockSize);
        }
        System.out.println(debugString + "  file position: " + currentFilePosition + " nSamples: " + nSamples);
        return true;
    }
}
