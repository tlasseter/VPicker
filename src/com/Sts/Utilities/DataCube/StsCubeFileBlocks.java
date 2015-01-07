package com.Sts.Utilities.DataCube;

import com.Sts.DBTypes.*;
import com.Sts.IO.StsAbstractFile;
import com.Sts.MVC.*;
import com.Sts.Reflect.*;
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

/** Defines a set of blocks for a cube of data.  These blocks could be in row, column, or slice order as defined by dir.
 *  For a given direction, though, think of it as row-ordered relative to that direction, so that data increases in order
 *  from slice to column to row, meaning that slice index increases fastest, then column, then row.  So if the dir is column or XDIR,
 *  then a column is a row and a row is a column in this class.  If the dir is slice or ZDIR, then a slice is a row and
 *  a row is a slice.
 */
public class StsCubeFileBlocks extends StsFileBlocks
{
	int dir;
	int nSlices;
	int nRows;
	int nCols;
	long nTotalPlanes;
	long nSamplesPerPlane;
	long nBytesPerPlane;
	public long nPlanesPerBlock;
	long nBytesPerBlock;
	long nSamplesPerBlock;

	static final int XDIR = StsParameters.XDIR;
	static final int YDIR = StsParameters.YDIR;
	static final int ZDIR = StsParameters.ZDIR;

	static public final float nullValue = StsParameters.nullValue;
	static public final byte nullByte = StsParameters.nullByte;

	static final boolean debug = false;

    /** This constructor divides the cube file into blocks with the number of planes per block determined by the maxBlockSize
     * @param dir direction of planes (X, Y, or Z
     * @param nRows number of row elements (incrementing Y)
     * @param nCols number of col elements (incrementing X)
     * @param nSlices number of slice elements (ncrementing Z)
     * @param seismicBoundingBox seismicBox which defines directory where the file to be used or constructed exists
     * @param filename name of the file
     * @param mode Read "r" or read-write "rw"
     * @param nBytesPerSample 1 for bytes and 4 for floats
     * @param blocksMemoryManager memory manager handling acquiring and releasing blocks
     */
    public StsCubeFileBlocks(int dir, int nRows, int nCols, int nSlices, StsSeismicBoundingBox seismicBoundingBox,
                             String filename, String mode, int nBytesPerSample, StsBlocksMemoryManager blocksMemoryManager)
                             throws IOException
	{
		this.dir = dir;
		this.nRows = nRows;
		this.nCols = nCols;
		this.nSlices = nSlices;
        fileSize = (long)nRows*nCols*nSlices*nBytesPerSample;
        initialize(seismicBoundingBox, filename, mode, nBytesPerSample, blocksMemoryManager);

		if(dir == YDIR)
		{
			nTotalPlanes = nRows;
			nSamplesPerPlane = nSlices * nCols;
		}
		else if(dir == XDIR)
		{
			nTotalPlanes = nCols;
			nSamplesPerPlane = nSlices * nRows;
		}
		else // dirNo == ZDIR
		{
			nTotalPlanes = nSlices;
			nSamplesPerPlane = nRows * nCols;
		}

	    long nTotalBytes = (long)nRows * (long)nCols * (long)nSlices * nBytesPerSample;

		nBytesPerPlane = nSamplesPerPlane * nBytesPerSample;

		if(nTotalBytes < maxBlockSize)
		{
			nPlanesPerBlock = nTotalPlanes;
			nTotalBlocks = 1;
			blocks = new StsFileBlock[1];
			blocks[0] = new StsFileBlock(this, 0, 0, nTotalBytes);
			nBytesPerBlock = nTotalBytes;
		}
		else
		{
            nPlanesPerBlock = maxBlockSize/nBytesPerPlane;
//            nPlanesPerBlock = Math.round(approxBlockSize / (nBytesPerPlane));
			nPlanesPerBlock = Math.max(1, nPlanesPerBlock);
			nTotalBlocks = StsMath.ceiling((float)nTotalPlanes / nPlanesPerBlock);
			nSamplesPerBlock = nPlanesPerBlock * nSamplesPerPlane;
			nBytesPerBlock = nSamplesPerBlock * nBytesPerSample;
			blocks = new StsFileBlock[nTotalBlocks];
			long offset = 0;
			for(int n = 0; n < nTotalBlocks-1; n++)
			{
				blocks[n] = new StsFileBlock(this, n, offset, nBytesPerBlock);
				offset += nBytesPerBlock;
			}
			long nBytesRemaining = nTotalBytes - (nTotalBlocks-1)*nBytesPerBlock;
			if(nBytesRemaining > 0)
			{
				blocks[nTotalBlocks-1] = new StsFileBlock(this, nTotalBlocks-1, offset, nBytesRemaining);
			}
		}
	}
    /** This constructor honors an input nPlanesPerBlock (often 1), so the allocated block is exactly a plane
     *  We may also wish this fileBlock set matches an already constructed fileBlocks set
     *  such as when we have defined a float set and want a byte set to have a matching planes set when they are used simultaneously in
     *  computing values
     * @param dir direction of planes (X, Y, or Z
     * @param nRows number of row elements (incrementing Y)
     * @param nCols number of col elements (incrementing X)
     * @param nSlices number of slice elements (ncrementing Z)
     * @param filename name of the file
     * @param mode Read "r" or read-write "rw"
     * @param nBytesPerSample 1 for bytes and 4 for floats
     * @param blocksMemoryManager memory manager handling acquiring and releasing blocks
     */
    public StsCubeFileBlocks(int dir, int nRows, int nCols, int nSlices, int nPlanesPerBlock, StsSeismicBoundingBox seismicBoundingBox,
                             String filename, String mode, int nBytesPerSample, StsBlocksMemoryManager blocksMemoryManager)
                             throws IOException
	{
		this.dir = dir;
		this.nRows = nRows;
		this.nCols = nCols;
		this.nSlices = nSlices;
        this.nPlanesPerBlock = nPlanesPerBlock;

        fileSize = (long)nRows*nCols*nSlices*nBytesPerSample;
        initialize(seismicBoundingBox, filename, mode, nBytesPerSample, blocksMemoryManager);

		if(dir == YDIR)
		{
			nTotalPlanes = nRows;
			nSamplesPerPlane = nSlices * nCols;
		}
		else if(dir == XDIR)
		{
			nTotalPlanes = nCols;
			nSamplesPerPlane = nSlices * nRows;
		}
		else // dirNo == ZDIR
		{
			nTotalPlanes = nSlices;
			nSamplesPerPlane = nRows * nCols;
		}

	    long nTotalBytes = (long)nRows * (long)nCols * (long)nSlices * nBytesPerSample;
		nBytesPerPlane = nSamplesPerPlane * nBytesPerSample;
        nTotalBlocks = StsMath.ceiling((float)nTotalPlanes / nPlanesPerBlock);
        nSamplesPerBlock = nPlanesPerBlock * nSamplesPerPlane;
        nBytesPerBlock = nSamplesPerBlock * nBytesPerSample;
        if(nBytesPerBlock > maxBlockSize)
        {
            new StsMessage(StsModel.getCurrentModel().win3d,  StsMessage.WARNING,
                    "Block allocation for file " + filename + " exceeds max allowed block size: " + maxBlockSize + "\n. Continuing...");
        }
        blocks = new StsFileBlock[nTotalBlocks];
        long offset = 0;
        for(int n = 0; n < nTotalBlocks-1; n++)
        {
            blocks[n] = new StsFileBlock(this, n, offset, nBytesPerBlock);
            offset += nBytesPerBlock;
        }
        long nBytesRemaining = nTotalBytes - (nTotalBlocks-1)*nBytesPerBlock;
        if(nBytesRemaining > 0)
        {
            blocks[nTotalBlocks-1] = new StsFileBlock(this, nTotalBlocks-1, offset, nBytesRemaining);
        }
	}
     /** This constructor honors an input nPlanesPerBlock (often 1), so the allocated block is exactly a plane
     *  We may also wish this fileBlock set matches an already constructed fileBlocks set
     *  such as when we have defined a float set and want a byte set to have a matching planes set when they are used simultaneously in
     *  computing values
     * @param dir direction of planes (X, Y, or Z
     * @param nRows number of row elements (incrementing Y)
     * @param nCols number of col elements (incrementing X)
     * @param nSlices number of slice elements (ncrementing Z)
     * @param file file reference from jar file
     * @param nBytesPerSample 1 for bytes and 4 for floats
     * @param blocksMemoryManager memory manager handling acquiring and releasing blocks
     */
    public StsCubeFileBlocks(int dir, int nRows, int nCols, int nSlices, int nPlanesPerBlock, File file,
                             int nBytesPerSample, StsBlocksMemoryManager blocksMemoryManager)
                             throws IOException
	{
		this.dir = dir;
		this.nRows = nRows;
		this.nCols = nCols;
		this.nSlices = nSlices;
        this.nPlanesPerBlock = nPlanesPerBlock;

        fileSize = (long)nRows*nCols*nSlices*nBytesPerSample;
        initialize(file, nBytesPerSample, blocksMemoryManager);

		if(dir == YDIR)
		{
			nTotalPlanes = nRows;
			nSamplesPerPlane = nSlices * nCols;
		}
		else if(dir == XDIR)
		{
			nTotalPlanes = nCols;
			nSamplesPerPlane = nSlices * nRows;
		}
		else // dirNo == ZDIR
		{
			nTotalPlanes = nSlices;
			nSamplesPerPlane = nRows * nCols;
		}

	    long nTotalBytes = (long)nRows * (long)nCols * (long)nSlices * nBytesPerSample;
		nBytesPerPlane = nSamplesPerPlane * nBytesPerSample;
        nTotalBlocks = StsMath.ceiling((float)nTotalPlanes / nPlanesPerBlock);
        nSamplesPerBlock = nPlanesPerBlock * nSamplesPerPlane;
        nBytesPerBlock = nSamplesPerBlock * nBytesPerSample;
        if(nBytesPerBlock > maxBlockSize)
        {
            new StsMessage(StsModel.getCurrentModel().win3d,  StsMessage.WARNING,
                    "Block allocation for file " + file.getName() + " exceeds max allowed block size: " + maxBlockSize + "\n. Continuing...");
        }
        blocks = new StsFileBlock[nTotalBlocks];
        long offset = 0;
        for(int n = 0; n < nTotalBlocks-1; n++)
        {
            blocks[n] = new StsFileBlock(this, n, offset, nBytesPerBlock);
            offset += nBytesPerBlock;
        }
        long nBytesRemaining = nTotalBytes - (nTotalBlocks-1)*nBytesPerBlock;
        if(nBytesRemaining > 0)
        {
            blocks[nTotalBlocks-1] = new StsFileBlock(this, nTotalBlocks-1, offset, nBytesRemaining);
        }
	}

    public boolean getRowCubeFloatValues(int row, int col, int sliceMin, int sliceMax, float[] values)
	{
		int position = getRowCubePosition(row, col, sliceMin, sliceMax);
		if(position == -1)
		{
			return false;
		}
		int nValues = sliceMax - sliceMin + 1;
		int nBlock = (int)(row / nPlanesPerBlock);
		MappedByteBuffer blockBuffer = getBlockBuffer(nBlock, FileChannel.MapMode.READ_ONLY);
		if(blockBuffer == null)
		{
			return false;
		}
		return getFloatValues(blockBuffer, position, nValues, values);
	}

	private MappedByteBuffer getReadOnlyBlock(int nBlock)
	{
		return getBlockBuffer(nBlock, FileChannel.MapMode.READ_ONLY);
	}

	private MappedByteBuffer getReadWriteBlock(int nBlock)
	{
		return getBlockBuffer(nBlock, FileChannel.MapMode.READ_WRITE);
	}

	int getRowCubePosition(int row, int col, int sliceMin, int sliceMax)
	{
		if(row < 0 || row >= nRows)
		{
			return -1;
		}
		if(col < 0 || col >= nCols)
		{
			return -1;
		}

		int nValues = sliceMax - sliceMin + 1;
		if(nValues <= 0)
		{
			return -1;
		}

		int nBlockRow = (int)(row % nPlanesPerBlock);
		int index = col * nSlices + sliceMin;
		return (int)(nBlockRow * nSamplesPerPlane + index);
	}

	int getRowCubePosition(int row, int col, int slice)
	{
        if(!insideVolume(row, col, slice)) return -1;
		int nBlockPlane = (int)(row % nPlanesPerBlock);
		int planeIndex = col * nSlices + slice;
		return (int)(nBlockPlane * nSamplesPerPlane + planeIndex);
	}

    int getColCubePosition(int row, int col, int slice)
    {
        if(!insideVolume(row, col, slice)) return -1;
        int nBlockPlane = (int)(col % nPlanesPerBlock);
        int planeIndex = row * nSlices + slice;
        return (int)(nBlockPlane * nSamplesPerPlane + planeIndex);
    }

    int getSliceCubePosition(int row, int col, int slice)
    {
        if(!insideVolume(row, col, slice)) return -1;
        int nBlockPlane = (int)(slice % nPlanesPerBlock);
        int planeIndex = row * nCols + col;
        return (int)(nBlockPlane * nSamplesPerPlane + planeIndex);
    }

    public boolean insideVolume(int row, int col, int slice)
    {
        if(row < 0 || row >= nRows) return false;
        if(col < 0 || col >= nCols) return false;
        if(slice < 0 || slice >= nSlices) return false;
        return true;
    }

    public float getRowCubeFloatValue(int row, int col, int slice)
	{
		int position = getRowCubePosition(row, col, slice);
		if(position == -1)
		{
			return nullValue;
		}
		long nBlock = row / nPlanesPerBlock;
		MappedByteBuffer blockBuffer = getBlockBuffer(nBlock, FileChannel.MapMode.READ_ONLY);
		if(blockBuffer == null)
		{
			return nullValue;
		}
		return getFloatValue(blockBuffer, position);
	}

	public float getRowCubeFloatValueDebug(int row, int col, int slice)
	{
		int position = getRowCubePosition(row, col, slice);
		if(position == -1)
		{
			return nullValue;
		}
		long nBlock = row / nPlanesPerBlock;
		MappedByteBuffer blockBuffer = getBlockBuffer(nBlock, FileChannel.MapMode.READ_ONLY);
		if(blockBuffer == null)
		{
			return nullValue;
		}
		long filePosition = nBlock * nSamplesPerBlock + position;
		long correctPosition = nSamplesPerPlane * row + nSlices * col +
			slice;
		System.out.println("    file position " + filePosition +
						   " correctPosition " + correctPosition);
		return getFloatValue(blockBuffer, position);
	}

	public boolean getRowCubeByteValues(int row, int col, int sliceMin, int sliceMax, byte[] values)
	{
		int position = getRowCubePosition(row, col, sliceMin, sliceMax);
		if(position == -1)
		{
			return false;
		}
		int nValues = sliceMax - sliceMin + 1;
		long nBlock = row / nPlanesPerBlock;
		MappedByteBuffer blockBuffer = getBlockBuffer(nBlock, FileChannel.MapMode.READ_ONLY);
		if(blockBuffer == null)
		{
			return false;
		}
		return getByteValues(blockBuffer, position, nValues, values);
	}

	public boolean getRowCubePlaneValues(int row, byte[] values)
	{
		int position = getRowCubePosition(row, 0, 0);
		if(position == -1)return false;
		long nBlock = row / nPlanesPerBlock;
		MappedByteBuffer blockBuffer = getBlockBuffer(nBlock, FileChannel.MapMode.READ_ONLY);
		if(blockBuffer == null)return false;
		int nValues = nCols * nSlices;
		return getByteValues(blockBuffer, position, nValues, values);
	}

	public boolean getRowCubePlaneValues(int row, int[] values)
	{
		int position = getRowCubePosition(row, 0, 0);
		if(position == -1)return false;
		long nBlock = row / nPlanesPerBlock;
		MappedByteBuffer blockBuffer = getBlockBuffer(nBlock, FileChannel.MapMode.READ_ONLY);
		if(blockBuffer == null)return false;
		int nValues = nCols * nSlices;
		return getUnsignedIntValues(blockBuffer, position, nValues, values);
	}

	public boolean putRowCubePlaneValues(int row, int nInputValues, int inputOffset, byte[] inputValues)
	{
        if(dir != YDIR)
        {
            StsException.systemError("Cannot put rowCube inputValues as this is not a row cube. Cube type is " + StsParameters.coorLabels[dir]);
            return false;
        }
        try
        {
            int rowBlockOutputOffset = getRowCubePosition(row, 0, 0);
            if(rowBlockOutputOffset == -1)return false;
            long nBlock = row / nPlanesPerBlock;
            MappedByteBuffer blockBuffer = getBlockBuffer(nBlock, FileChannel.MapMode.READ_WRITE);
            if(blockBuffer == null)return false;
            putByteValues(blockBuffer, rowBlockOutputOffset, nInputValues, inputOffset, inputValues);
            return true;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "putRowCubePlaneValues(row, byte[])", e);
            return false;
        }
    }

    /** Puts byte inputValues into col cube mappedByteBuffer beginning at a trace locted at specified row and column */
    public boolean putColCubePlaneValues(int row, int col, int nInputValues, int inputOffset, byte[] inputValues)
    {
        if(dir != XDIR)
        {
            StsException.systemError("Cannot put colCube inputValues as this is not a col cube. Cube type is " + StsParameters.coorLabels[dir]);
            return false;
        }
        try
        {
            int colBlockOutputOffset = getColCubePosition(row, col, 0);
            if(colBlockOutputOffset == -1)return false;
            int nBlock = (int)(col / nPlanesPerBlock);
            MappedByteBuffer blockBuffer = getBlockBuffer(nBlock, FileChannel.MapMode.READ_WRITE);
            if(blockBuffer == null)return false;
            putByteValues(blockBuffer, colBlockOutputOffset, nInputValues, inputOffset, inputValues);
            // we may not be thru with this block, but unlock it anyways in case we need the memory
            unlockBlock(nBlock);
            return true;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "putColCubePlaneValues(row, byte[])", e);
            return false;
        }
    }

    public boolean putSliceCubePlaneValues(int row, int slice, int nInputValues, int inputOffset, byte[] inputValues)
    {
        if(dir != ZDIR)
        {
            StsException.systemError("Cannot put colCube inputValues as this is not a col cube. Cube type is " + StsParameters.coorLabels[dir]);
            return false;
        }
        try
        {
            int sliceBlockOutputOffset = getSliceCubePosition(row, 0, slice);
            if(sliceBlockOutputOffset == -1)return false;
            long nBlock = slice / nPlanesPerBlock;
            MappedByteBuffer blockBuffer = getBlockBuffer(nBlock, FileChannel.MapMode.READ_WRITE);
            if(blockBuffer == null)return false;
            putByteValues(blockBuffer, sliceBlockOutputOffset, nInputValues, inputOffset, inputValues);
            return true;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "putSliceCubePlaneValues(row, byte[])", e);
            return false;
        }
    }
	/*
	 boolean getValues(int row, int col, int sliceMin, int sliceMax, byte[] values)
	  {
	   int length = sliceMax - sliceMin + 1;
	   if(length <= 0) return false;

	   byte[] planeData;

	   int n = 0;
	   int index = row*nCols + col;
	   for(int slice = sliceMin; slice <= sliceMax; slice++)
	 values[n++] = readValue(slice, index);
	   return true;
	  }
	 */

	public float readFloatValue(int row, int col, int slice)
	{
		int planeIndex = -1;
		long nBlock = -1;
		long nBlockPlane = -1;
		try
		{
			if(dir == XDIR)
			{
				nBlock = col / nPlanesPerBlock;
				nBlockPlane = col % nPlanesPerBlock;
				planeIndex = row * nSlices + slice;
			}
			if(dir == YDIR)
			{
				nBlock = row / nPlanesPerBlock;
				nBlockPlane = row % nPlanesPerBlock;
				planeIndex = col * nSlices + slice;
			}
			else // dirNo == ZDIR
			{
				nBlock = slice / nPlanesPerBlock;
				nBlockPlane = slice % nPlanesPerBlock;
				planeIndex = row * nCols + col;
			}
			MappedByteBuffer blockBuffer = getBlockBuffer(nBlock, FileChannel.MapMode.READ_ONLY);
			if(blockBuffer == null)
			{
				return nullValue;
			}
			FloatBuffer floatBuffer = blockBuffer.asFloatBuffer();
			int nValue = 0;
            int position = (int)(nBlockPlane * nSamplesPerPlane + planeIndex);
            return floatBuffer.get(position);
		}
		catch(Exception e)
		{
			StsException.outputException(
				"StsSeismicVolume.readValue() failed. Couldn't read value for row " +
				row + " col " + col + "slice " + slice + " index: " + planeIndex + " nBlock: " +
				nBlock, e, StsException.WARNING);
			return nullValue;
		}
	}

	public float[] readFloatPlane(int nPlane)
	{
		long nBlock = -1;
		try
		{
			nBlock = nPlane / nPlanesPerBlock;
			long nBlockPlane = nPlane % nPlanesPerBlock;
			MappedByteBuffer blockBuffer = getBlockBuffer(nBlock, FileChannel.MapMode.READ_ONLY);
			if(blockBuffer == null)return null;
			FloatBuffer floatBuffer = blockBuffer.asFloatBuffer();
            int position = (int)(nBlockPlane * nSamplesPerPlane);
            floatBuffer.position(position);
			float[] floats = new float[(int)nSamplesPerPlane];
            floatBuffer.get(floats);
			return floats;
		}
		catch(Exception e)
		{
			StsException.outputException(
				"StsSeismicVolume.readValue() failed. Couldn't read plane " + nPlane + " from volume " + pathname +
				" block " + nBlock, e, StsException.WARNING);
			return null;
		}
	}

    public boolean readFloatPlane(int nPlane, float[] floats)
	{
		long nBlock = -1;
		try
		{
			nBlock = nPlane / nPlanesPerBlock;
			long nBlockPlane = nPlane % nPlanesPerBlock;
			MappedByteBuffer blockBuffer = getBlockBuffer(nBlock, FileChannel.MapMode.READ_ONLY);
			if(blockBuffer == null)return false;
			FloatBuffer floatBuffer = blockBuffer.asFloatBuffer();
            int position = (int)(nBlockPlane * nSamplesPerPlane);
            floatBuffer.position(position);
            if(floats == null)
            {
                StsException.systemError(this, "readFloatPlane", "floats array is null. Cannot read data into it.");
                return false;
            }
            else if(floats.length < nSamplesPerPlane)
            {
                StsException.systemError(this, "readFloatPlane", "floats array is short. Length is " + floats.length + ", should be: " + nSamplesPerPlane);
                return false;
            }
            floatBuffer.get(floats, 0, (int)nSamplesPerPlane);
			return true;
		}
		catch(Exception e)
		{
			StsException.outputWarningException(this, "readFloatPlane", "Couldn't read plane " + nPlane + " from volume " + pathname + " block " + nBlock, e);
			return false;
		}
	}

    public boolean readBlockFloats(int nBlock, float[] floats)
	{
		try
		{
			MappedByteBuffer blockBuffer = getBlockBuffer(nBlock, FileChannel.MapMode.READ_ONLY);
			if(blockBuffer == null)return false;
			FloatBuffer floatBuffer = blockBuffer.asFloatBuffer();
            floatBuffer.get(floats);
            return true;
		}
		catch(Exception e)
		{
			StsException.outputWarningException(this, "readBlockFloats", "Couldn't read block " + nBlock + " from volume " + pathname, e);
			return false;
		}
	}

    public boolean writeBlockFloats(int nBlock, float[] floats)
	{
		try
		{
			MappedByteBuffer blockBuffer = getBlockBuffer(nBlock, FileChannel.MapMode.READ_WRITE);
			if(blockBuffer == null)return false;
			FloatBuffer floatBuffer = blockBuffer.asFloatBuffer();
            floatBuffer.put(floats);
            return true;
		}
		catch(Exception e)
		{
			StsException.outputWarningException(this, "writeBlockFloats", "Couldn't read block " + nBlock + " from volume " + pathname, e);
			return false;
		}
	}

    public byte readByteValue(int row, int col, int slice)
	{
		int planeIndex = -1;
		long nBlock = -1;
		long nBlockPlane = -1;
		try
		{
			if(dir == XDIR)
			{
				nBlock = col / nPlanesPerBlock;
				nBlockPlane = col % nPlanesPerBlock;
				planeIndex = row * nSlices + slice;
			}
			if(dir == YDIR)
			{
				nBlock = row / nPlanesPerBlock;
				nBlockPlane = row % nPlanesPerBlock;
				planeIndex = col * nSlices + slice;
			}
			else // dirNo == ZDIR
			{
				nBlock = slice / nPlanesPerBlock;
				nBlockPlane = slice % nPlanesPerBlock;
				planeIndex = row * nCols + col;
			}
			MappedByteBuffer blockBuffer = getBlockBuffer(nBlock, FileChannel.MapMode.READ_ONLY);
			if(blockBuffer == null)
			{
				return nullByte;
			}
            int position = (int)(nBlockPlane * nSamplesPerPlane + planeIndex);
            return blockBuffer.get(position);
		}
		catch(Exception e)
		{
			StsException.outputException(
				"StsSeismicVolume.readValue() failed. Couldn't read value for row " +
				row + " col " + col + "slice " + slice + " index: " + planeIndex + " nBlock: " +
				nBlock, e, StsException.WARNING);
			return nullByte;
		}
	}

	public byte[] readBytePlane(int nPlane)
	{
		long nBlock = -1;
		try
		{
			nBlock = nPlane / nPlanesPerBlock;
			long nBlockPlane = nPlane % nPlanesPerBlock;
			MappedByteBuffer blockBuffer = getBlockBuffer(nBlock, FileChannel.MapMode.READ_ONLY);
			if(blockBuffer == null)return null;
            int position = (int)(nBlockPlane * nSamplesPerPlane);
            blockBuffer.position(position);
			byte[] bytes = new byte[(int)nSamplesPerPlane];
			blockBuffer.get(bytes);
			return bytes;
		}
		catch(Exception e)
		{
			StsException.outputException(
				"StsSeismicVolume.readBytePlane() failed. Couldn't read plane " + nPlane + " from volume " + pathname +
				" block " + nBlock, e, StsException.WARNING);
			return null;
		}
	}

	public ByteBuffer getByteBufferPlane(int nPlane, FileChannel.MapMode mapMode)
	{
		long nBlock = -1;
		try
		{
			nBlock = nPlane / nPlanesPerBlock;
			long nBlockPlane = nPlane % nPlanesPerBlock;
			MappedByteBuffer blockBuffer = getBlockBuffer(nBlock, mapMode);
			if(blockBuffer == null)return null;
            int position = (int)(nBlockPlane * nBytesPerPlane);
            blockBuffer.position(position);
			ByteBuffer planeBuffer = blockBuffer.slice();
			return planeBuffer;
//			return blockBuffer;
		}
		catch(Exception e)
		{
			StsException.outputWarningException(this, "readByteBufferPlane",
                    "Couldn't read plane " + nPlane + " from volume " + pathname + " block " + nBlock, e);
			return null;
		}
	}

	public void unlockPlane(int nPlane)
	{
		int nBlock = (int)(nPlane/nPlanesPerBlock);
        unlockBlock(nBlock);
    }

    public Object calculateDataPlane(int nPlane_, StsMethod calculateDataPlaneMethod_)
	{
		// don't calculate while button is down or being dragged
//            boolean isDragging = currentModel.win3d.getCursor3d().getIsDragging();
//            System.out.println("line2d.calculateDataPlane() dirNo: " + dirNo + " plane: " + nPlane_ + " isDragging: " + isDragging);
		if(StsObject.getCurrentModel().getCursor3d().getIsDragging(dir))
		{
			return null;
		}

		final int nPlane = nPlane_;
		final byte[] planeData = new byte[(int)nBytesPerPlane];
		final StsMethod calculateDataPlaneMethod = calculateDataPlaneMethod_;

		/* ideally after calculating we would store plane in temp file; for now always recompute
			if (planes[nPlane] != null)
			{
		 return planes[nPlane];
			}
		 */
		Runnable calcRunnable = new Runnable()
		{
			public void run()
			{
				if(calculateDataPlane(nPlane, planeData, calculateDataPlaneMethod))
				{
					StsModel model = StsObject.getCurrentModel();
					model.getCursor3d().clearTextureDisplays(dir);
					model.win3dDisplay();
				}
			}
		};

		Thread calcThread = new Thread(calcRunnable);
		calcThread.start();

		return planeData;
	}

	public boolean calculateDataPlane(int nPlane, byte[] planeData, StsMethod calculateDataPlaneMethod)
	{
		try
		{
//				setPlaneData(nPlane, planeData);
			if(debug)
			{
				System.out.println("calculateDataPlane called with planeData[0] " + planeData[0] + " dirNo " + dir + " nPlane " + nPlane);
			}
			calculateDataPlaneMethod.invokeInstanceMethod(new Object[] { planeData, new Integer(dir), new Integer(nPlane) });
//				setPlaneData(nPlane, planeData);
			/*
			 if (file == null) return false;
			 file.seek( (long) nPlanePoints * nPlane);
			 file.write(planeData);return
			 return
			 */
			return true;
		}

		catch(Exception e)
		{
			StsException.outputException("StsSeismicVolume.calculateDataPlane() failed.", e, StsException.WARNING);
			return false;
		}
	}
/*
	public void deleteAllBlocks()
	{
		System.out.println("StsSeismicVoume.deleteAllBlocks() freeing " + nTotalBlocks + " of virtual memory rowBlocks.");
		for(int n = 0; n < nTotalBlocks; n++)
		{
			StsFileBlock block = blocks[n];
			block.clear();
		}
		memoryUsed = 0;
	}
*/
	public boolean getSliceBlockTraceValues(int row, int col, int sliceMin, int sliceMax, byte[] values)
	{
		int nSliceValues = sliceMax - sliceMin + 1;
		if(nSliceValues <= 0)
		{
			return false;
		}

		int index = row * nCols + col;
		long nBlock = sliceMin / nPlanesPerBlock;
		MappedByteBuffer blockByteBuffer = this.getBlockBuffer(nBlock, FileChannel.MapMode.READ_ONLY);
		if(blockByteBuffer == null)
		{
			return false;
		}
		int nValue = 0;
		long nBlockSlice = sliceMin % nPlanesPerBlock;
		long location = nBlockSlice*nSamplesPerPlane + index;
		for(int nSlice = sliceMin; nSlice <= sliceMax; nSlice++, nValue++, nBlockSlice++, location += nSamplesPerPlane)
		{
			if(nBlockSlice >= nPlanesPerBlock)
			{
				nBlock++;
				blockByteBuffer = getBlockBuffer(nBlock, FileChannel.MapMode.READ_ONLY);
				if(blockByteBuffer == null)
				{
					return false;
				}
				nBlockSlice = 0;
				location = index;
			}
			values[nValue] = blockByteBuffer.get((int)location);
		}
		return true;
	}

	public boolean getRowBlockTraceValues(int row, int col, int sliceMin, int sliceMax, byte[] values)
	{
		int nSliceValues = sliceMax - sliceMin + 1;
		if(nSliceValues <= 0) return false;
		long nBlock = row / nPlanesPerBlock;
		long nBlockRow = row % nPlanesPerBlock;
		try
		{
			MappedByteBuffer blockByteBuffer = getBlockBuffer(nBlock, FileChannel.MapMode.READ_ONLY);
			if(blockByteBuffer == null)return false;
			long location = nBlockRow * nSamplesPerPlane + col * nSlices + sliceMin;
			blockByteBuffer.position((int)location);
			blockByteBuffer.get(values);
			return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsCubeFileBlocks.getRowBlockValues() failed.", e, StsException.WARNING);
			return false;
		}
	}

	public boolean getPlane(int nPlane, byte[] values)
	{
		long nBlock, nBlockPlane;

		try
		{
			nBlock = nPlane / nPlanesPerBlock;
			MappedByteBuffer blockByteBuffer = getBlockBuffer(nBlock, FileChannel.MapMode.READ_ONLY);
			if(blockByteBuffer == null)
			{
				return false;
			}
			nBlockPlane = nPlane % nPlanesPerBlock;
            int position = (int)(nBlockPlane * nSamplesPerPlane);
            blockByteBuffer.position(position);
            blockByteBuffer.get(values, 0, (int)nSamplesPerPlane);
			return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsSeismicVolume.FileMapSliceBlocks.getSlice() failed.", e, StsException.WARNING);
			return false;
		}
	}

	public int getPlanePosition(int nPlane)
	{
		return (int)((nPlane % nPlanesPerBlock) * nSamplesPerPlane);
	}
/*
	public MappedByteBuffer getPlaneBuffer(int nPlane)
	{
		int nBlock, nBlockPlane;

		try
		{
			nBlock = nPlane / nPlanesPerBlock;
			MappedByteBuffer blockByteBuffer = getReadOnlyBlock(nBlock);
			if(blockByteBuffer == null)
			{
				return null;
			}
			nBlockPlane = nPlane % nPlanesPerBlock;
			blockByteBuffer.position(nBlockPlane * nSamplesPerPlane);
			return blockByteBuffer;
		}
		catch(Exception e)
		{
			StsException.outputException("StsSeismicVolume.FileMapSliceBlocks.getSlice() failed.", e, StsException.WARNING);
			return null;
		}
	}
*/
	public void printMemorySummary()
	{
		int nBlocksInMemory = 0;
		for(int n = 0; n < nTotalBlocks; n++)
			if(blocks[n].byteBuffer != null) nBlocksInMemory++;
		StsMessageFiles.infoMessage("Blocks in memory: " + nBlocksInMemory + ". Memory used: " + memoryUsed);
	}
	public static void main(String[] args)
	{
		try
		{
			String directory = System.getProperty("user.dirNo");
			String filename = "testWriteBlocks";
			StsBlocksMemoryManager memoryManager = new StsBlocksMemoryManager(20);
            StsSeismicBoundingBox seismicBox = new StsSeismicBoundingBox(0, 100, 0, 50, 0, 25);
            StsCubeFileBlocks cubeFileByteBlocks = new StsCubeFileBlocks(0, 400, 400, 2000, seismicBox, filename, "rw", 1, memoryManager);
			ByteBuffer stackedDataTest = writeTestByteData(cubeFileByteBlocks, 0, 100, 400, 400, 2000);
			stackedDataTest.rewind();
			readTestByteData(cubeFileByteBlocks, 0, 100, 400, 400, 2000);
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}

	static public ByteBuffer writeTestByteData(StsCubeFileBlocks fileByteBlocks, int dir, int nPlane, int nRows, int nCols, int nSlices)
	{
		ByteBuffer byteBuffer = null;
		int row = -1, col = -1;
		try
		{
			byteBuffer = fileByteBlocks.getByteBufferPlane(nPlane, FileChannel.MapMode.READ_WRITE);
			System.out.println("  ByteBuffer capacity: " + byteBuffer.capacity());
			byte[] data = new byte[nSlices];
			{
				//            progressBarDialog.setLabelText("Stacking Crossline #" + this.getColNumFromCol(index));
				//            progressBarDialog.setProgressMax(nRows);
				//            progressBarDialog.pack();
				//			  progressBarDialog.setVisible(true);
				int n = 0;
				col = nPlane;
				for(row = 0; row < nRows; row++, n++)
				{
					data[0] = (byte)row;
					System.out.println("ByteBuffer.put data trace for row " + row + " col " + col + " data[0] " + data[0] + "  ByteBuffer position: " + byteBuffer.position() + " isDirect " + byteBuffer.isDirect());
					byteBuffer.put(data);
				}
			}
			return byteBuffer;
		}
		catch(Exception e)
		{
			if(byteBuffer != null)
				StsException.systemError("StsPreStackLineSet3d.getStackedByteData() failed. \n" +
					"byteBuffer capacity: " + byteBuffer.capacity() + " position: " + byteBuffer.position());
			StsException.outputException("StsPreStackLineSet3d.getStackedByteData() failed.", e, StsException.WARNING);
			return null;
		}
	}
	static public void readTestByteData(StsCubeFileBlocks fileByteBlocks, int dir, int nPlane, int nRows, int nCols, int nSlices)
	{
		ByteBuffer byteBuffer = null;
		int row = -1, col = -1;
		try
		{
			byteBuffer = fileByteBlocks.getByteBufferPlane(nPlane, FileChannel.MapMode.READ_WRITE);
			System.out.println("  ByteBuffer capacity: " + byteBuffer.capacity());
			byte[] data = new byte[nSlices];
			{
				//            progressBarDialog.setLabelText("Stacking Crossline #" + this.getColNumFromCol(index));
				//            progressBarDialog.setProgressMax(nRows);
				//            progressBarDialog.pack();
				//			  progressBarDialog.setVisible(true);
				int n = 0;
				col = nPlane;
				for(row = 0; row < nRows; row++, n++)
				{
					byteBuffer.get(data);
					System.out.println("ByteBuffer.get data trace for row " + row + " col " + col + " data[0] " + data[0] + "  ByteBuffer position: " + byteBuffer.position() + " isDirect " + byteBuffer.isDirect());
				}
			}
		}
		catch(Exception e)
		{
			if(byteBuffer != null)
				StsException.systemError("StsPreStackLineSet3d.getStackedByteData() failed. \n" +
					"byteBuffer capacity: " + byteBuffer.capacity() + " position: " + byteBuffer.position());
			StsException.outputException("StsPreStackLineSet3d.getStackedByteData() failed.", e, StsException.WARNING);
		}
	}
/*
	public static void main(String[] args)
	{
		try
		{
			String pathname = System.getProperty("user.dirNo");
			StsFileChooser chooser = new StsFileChooser(null,
				"select a seismic volume", pathname, "seis3d.txt", true, false);
			chooser.show();
			File[] files = chooser.getFiles();
			if (files == null || files[0] == null)
			{
				System.out.println("No file selected.");
				System.exit(0);
			}
			File file = files[0];
			pathname = file.getPath();
			StsSeismicVolume volume = new StsSeismicVolume();
			StsParameterFile.readObjectFields(pathname, volume,
											  StsSeismicBoundingBox.class,
											  StsBoundingBox.class);
			testSliceBlocks(volume);
		}
		catch (Exception e)
		{
			StsException.outputException("StsSeismicVolume.main() failed.", e,
										 StsException.WARNING);
		}
		finally
		{
			System.exit(0);
		}
	}

	static void testSliceBlocks(StsSeismicVolume volume)
	{
		try
		{
			volume.allocateVolumes("r");
			StsCubeFileBlocks filesMapBlocks[] = volume.getFilesMapBlocks();
			for(int i = 0; i < 3; i++)
			{
				StsCubeFileBlocks fileMapBlock = filesMapBlocks[i];
				int nTotalBlocks = fileMapBlock.nTotalBlocks;
				int nPlanesPerBlock = fileMapBlock.nPlanesPerBlock;
				int nSamplesPerPlane = fileMapBlock.nSamplesPerPlane;
				byte[] values = new byte[nSamplesPerPlane];
				int nPlane = 0;
				for(int n = 0; n < nTotalBlocks; n++)
				{
					fileMapBlock.getPlane(nPlane, values);
					System.out.println("Read first value from plane " + nPlane + " in block " + n + " of " + nTotalBlocks + " blocks: " + values[0]);
					nPlane += nPlanesPerBlock;
				}
			}
		}
		catch (Exception e)
		{
			System.out.println(
				"Exception thrown by StsSeismicVolume.testSliceBlocks()");
			e.printStackTrace();
		}
	}
*/
}
