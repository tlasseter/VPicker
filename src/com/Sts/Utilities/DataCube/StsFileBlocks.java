package com.Sts.Utilities.DataCube;

import com.Sts.IO.StsAbstractFile;
import com.Sts.IO.StsJarEntryFile;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.nio.*;
import java.nio.channels.*;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */
abstract public class StsFileBlocks
{
	public int nTotalBlocks;
	int nBytesPerSample = 1;
	FileChannel channel = null;
    RandomAccessFile file = null;
	public StsFileBlock[] blocks = null;
	public String pathname;
    public String mode;
    long memoryUsed = 0;
    long fileSize;

    /** rough size of a memory-mapped block of seismic planes */

//	static public long approxBlockSize = StsMemoryAllocation.maxMapBufferSize;
	static public long maxBlockSize = StsMemoryAllocation.maxMapBufferSize;

     // small block testing
	//long approxBlockSize = 1000000;
	//long maxBlockSize = 2000000;

	static public StsBlocksMemoryManager blocksMemoryManager = null;

	static public final float nullValue = StsParameters.nullValue;
	static public final byte nullByte = StsParameters.nullByte;

	static final boolean debug = false;

	public StsFileBlocks()
	{
    }

	public void initialize(File file, int nBytesPerSample, StsBlocksMemoryManager blocksMemoryManager)
            throws FileNotFoundException, IOException
    {
		this.nBytesPerSample = nBytesPerSample;
		this.blocksMemoryManager = blocksMemoryManager;
		this.pathname = file.getAbsolutePath();
        this.mode = "r";
        initialize(file);
    }

	public void initialize(StsSeismicBoundingBox boundingBox, String filename, String mode, int nBytesPerSample, StsBlocksMemoryManager blocksMemoryManager)
            throws FileNotFoundException, IOException
    {
        if(mode.equals("r")) boundingBox.checkStsDirectoryForFilename(filename); 
		this.nBytesPerSample = nBytesPerSample;
		this.blocksMemoryManager = blocksMemoryManager;
        this.pathname = boundingBox.stsDirectory + filename;
        this.mode = mode;
        initialize();
    }

    /** Cases:
     *  1. file exists and readOnly: if right size, open and return; if wrong size complain and return
     *  2. file exists and readWrite: if right size, open and return; if wrong size, delete, create, and write 0 byte at end to claim space
     *  3. file doesn't exist and readOnly: error - throw and return
     *  4. file doesn't exist and readWrite: create, open and write 0 byte at end to claim space
     * @throws IOException
     */
    public void initialize() throws IOException
    {
        if(debug) System.out.println("Opening file and channel for " + pathname);
        File checkFile = new File(pathname);
        boolean fileExists = checkFile.exists();
        boolean readOnly =  mode.equals("r");
        if(fileExists)
        {
            long currentFileSize = checkFile.length();
            if(currentFileSize == fileSize)
            {
                openFile();
                return;
            }
            if(readOnly) // wrong size
            {
                /* For now, assume size doesn't matter and file is ok
                String message = "File " + pathname + " is wrong size.";
                StsMessage.printMessage(message);
                throw new FileNotFoundException(message);
                */
                openFile();
                return;
            }
            else // read-write file, wrong size, delete, create, and write zero byte to end 
            {
                checkFile.delete();
                openFile();
                writeEndByte(fileSize);
            }
        }
        else if(readOnly) // file doesn't exist and is read only: throw exception
        {
            String message = "Failed to find file " + pathname + ". Ignoring file and continuing.";
            StsMessage.printMessage(message);
            throw new FileNotFoundException(message);
        }
        else // file doesn't exist and is readWrite
        {
            openFile();
            writeEndByte(fileSize);
        }
    }

    /** Cases:
     *  1. file exists and readOnly: if right size, open and return; if wrong size complain and return
     *  2. file exists and readWrite: if right size, open and return; if wrong size, delete, create, and write 0 byte at end to claim space
     *  3. file doesn't exist and readOnly: error - throw and return
     *  4. file doesn't exist and readWrite: create, open and write 0 byte at end to claim space
     * @throws IOException
     */
    public void initialize(File file) throws IOException
    {
        if(debug) System.out.println("Opening file and channel for " + pathname);

        openFile(file);
        return;
    }

    private void openFile() throws IOException
    {
        file = new RandomAccessFile(pathname, mode);
        channel = file.getChannel();
    }

    private void openFile(File ifile) throws IOException
    {
        file = new RandomAccessFile(ifile, mode);
        channel = file.getChannel();
    }

    private void  writeEndByte(long fileSize)
    {
        try
        {
            file.seek(fileSize-1);
            file.write(0);
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "writeByte", "Failed to write zero byte at " + fileSize, e);
        }
    }


    // Added because the classInitialize that opens the channel never closes it and once we get to 1960 lines we are unable to load any more lines.
    // Moved the classInitialize to the getSuperGather method and close after successful read of block --- SAJ
    public void closeFile()
    {
        try
        {
            if(channel != null)
            {
                if(debug) System.out.println("Closing channel and file for " + pathname);
                channel.close();
                channel = null;
            }
            if(file != null)
            {
                file.close();
            }
        }
        catch(Exception e)
        {
            StsException.systemError("StsFileBlocks.closefile() failed for  " + pathname + ".");
        }
    }

    public FileChannel getChannel() { return channel; }

    public RandomAccessFile getFile() { return file; }


    public void unlockBlock(int nBlock)
    {
        if(nBlock < 0 || nBlock >= nTotalBlocks) return;
        blocks[nBlock].locked = false;
	}

    public void unlockAllBlocks()
    {
        for(int nBlock = 0; nBlock < nTotalBlocks; nBlock++)
            blocks[nBlock].locked = false;
	}

    public void unlockClearAllBlocks()
    {
        blocksMemoryManager.unlockAndClearBlocks(this);
	}

    static protected float getFloatValue(MappedByteBuffer blockBuffer, int position)
	{
		FloatBuffer floatBuffer = blockBuffer.asFloatBuffer();
		floatBuffer.position(position);
		return floatBuffer.get();
	}

	static public boolean getFloatValues(MappedByteBuffer blockBuffer, int position, int nValues, float[] values)
	{
		FloatBuffer floatBuffer = blockBuffer.asFloatBuffer();
		floatBuffer.position(position);
		floatBuffer.get(values, 0, nValues);
		return true;
	}

	static public boolean getByteValues(MappedByteBuffer blockBuffer, int position, int nValues, byte[] values)
	{
		blockBuffer.position(position);
		blockBuffer.get(values, 0, nValues);
		return true;
	}

	static public boolean getUnsignedIntValues(MappedByteBuffer blockBuffer, int position, int nValues, int[] values)
	{
		blockBuffer.position(position);
		byte[] bytes = new byte[nValues];
		blockBuffer.get(bytes);
		for(int i = 0; i < nValues; i++)
			values[i] = (int)(bytes[i] & 0xFF);
		return true;
	}

	static public FloatBuffer putFloatValue(MappedByteBuffer blockBuffer, int position, float value)
	{
		FloatBuffer floatBuffer = blockBuffer.asFloatBuffer();
		floatBuffer.position(position);
		return floatBuffer.put(value);
    }

	static public void putFloatValues(MappedByteBuffer blockBuffer, int position, int nValues, float[] values)
	{
		FloatBuffer floatBuffer = blockBuffer.asFloatBuffer();
		floatBuffer.position(position);
		floatBuffer.put(values, 0, nValues);
	}

	static public void putByteValues(MappedByteBuffer blockBuffer, int position, int nValues, int offset, byte[] values)
	{
		blockBuffer.position(position);
		blockBuffer.put(values, offset, nValues);
	}

	static public IntBuffer putUnsignedIntValues(MappedByteBuffer blockBuffer, int position, int[] values)
	{
        IntBuffer intBuffer = blockBuffer.asIntBuffer();
        intBuffer.position(position);
		return intBuffer.put(values);
	}

    public MappedByteBuffer getBlockBuffer(long nBlock, FileChannel.MapMode mapMode)
	{
		try
		{
            boolean readWrite = mapMode == FileChannel.MapMode.READ_WRITE;
            if(debug) System.out.println("Getting block " + nBlock + " readWrite: " + readWrite + " file: " + pathname);
            StsFileBlock block = blocks[(int)nBlock];
            if(block == null)
            {
			    StsException.systemError(this, "getBlock", "Block " + nBlock + " is null.");
                return null;
            }

            MappedByteBuffer blockBuffer = block.byteBuffer;

            
            if(blockBuffer != null)
            {
                // if we have this block and it is read-write or we want readOnly, rewind and return
                if(!blockBuffer.isReadOnly() || !readWrite)
                {
 //                    blockBuffer.rewind();
                     return blockBuffer;
                }
                // block is readOnly && we need readWrite: delete and remap it
                if(!blocksMemoryManager.unlockClearBlock(block))
                {
                    blockBuffer = null;
                    StsException.systemError(this, "getBlock", "Failed to clear read only block for read/write...continuing.");
                }
            }

			if(blockBuffer == null)
			{
				blocksMemoryManager.checkMemory();
				try
				{
                    if(channel == null) initialize();
                    blockBuffer = channel.map(mapMode, block.offset, block.size);
				}
				catch(IllegalArgumentException iae)
				{
					StsException.outputWarningException(this, "getBlock" + "Verify that you have enough disk space on " + this.pathname + "\n" +
                                            "Couldn't create block. IllegalArgumentException for block: " + nBlock, iae);
					return null;
				}
				catch(IOException ioe1)
				{
					boolean clearTheseBlocks = blocksMemoryManager.clearBlocks(this);
					boolean clearOtherBlocks = blocksMemoryManager.clearOtherBlocks(this);
					if(!clearTheseBlocks && !clearOtherBlocks)
					{
						StsException.systemError(this, "getBlock", "Could not free any other blocks.");
						return null;
					}

					try
					{
						blockBuffer = channel.map(mapMode, block.offset, block.size);
					}
					catch(IOException ioe2)
					{
						StsException.outputWarningException(this, "getBlock", "Could not create blockByteBuffer for block: " +
							nBlock + " offset: " + block.offset + " size: " + block.size + ". Error: "+ ioe2.getMessage(), ioe2);
						return null;
					}
				}
				block.byteBuffer = blockBuffer;
				if(debug)
				{
					System.out.println("mapped block " + pathname + "[" + nBlock + "]");
				}
				memoryUsed += block.size;
                block.locked = readWrite;
                blocksMemoryManager.addBlock(block);
			}
			else
			{
				blockBuffer.rewind();
			}
			return blockBuffer;
		}
		catch(Exception e)
		{
			StsException.outputWarningException(this, "getBlock", "Couldn't create block: " + nBlock, e);
			return null;
		}
	}

	static public void checkMemory(StsModel model)
	{
		StsBlocksMemoryManager blocksMemory = model.getProject().getBlocksMemoryManager();
		blocksMemory.checkMemory();
	}

	public boolean blockDeleted(long blockMemory)
	{
		memoryUsed -= blockMemory;
		return true;
	}

	public void deleteAllBlocks()
	{
		if(debug) StsException.systemDebug(this, "deleteAllBlocks", " freeing " + nTotalBlocks + " of virtual memory rowBlocks.");
		for(int n = 0; n < nTotalBlocks; n++)
		{
			StsFileBlock block = blocks[n];
			MappedByteBuffer byteBuffer = block.byteBuffer;
			if(byteBuffer != null) StsToolkit.clean(byteBuffer);
		}
		memoryUsed = 0;
	}

	public void printMemorySummary()
	{
		int nBlocksInMemory = 0;
		for(int n = 0; n < nTotalBlocks; n++)
			if(blocks[n].byteBuffer != null) nBlocksInMemory++;
		StsMessageFiles.infoMessage("Blocks in memory: " + nBlocksInMemory + ". Memory used: " + memoryUsed);
	}
}
