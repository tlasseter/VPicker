package com.Sts.Utilities.DataCube;

import com.Sts.Utilities.*;

import java.nio.*;

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
public class StsFileBlock
{
	public StsFileBlocks fileMapBlocks;
	public int nBlock;
	public long offset;
	public long size;
	public MappedByteBuffer byteBuffer = null;
	public boolean locked = false;

    public StsFileBlock(StsFileBlocks fileMapBlocks, int nBlock, long offset, long size)
	{
		this.fileMapBlocks = fileMapBlocks;
		this.nBlock = nBlock;
		this.offset = offset;
		this.size = size;
        if(StsFileBlocks.debug) StsException.systemDebug(this, "constructor", " for file: " + fileMapBlocks.pathname + " block: " + nBlock);
    }

    public boolean matches(StsFileBlocks fileMapBlocks, int nBlock)
    {
        return this.fileMapBlocks == fileMapBlocks && this.nBlock == nBlock;
    }

    public boolean unlockClear()
    {
        locked = false;
        return clear();
    }
    
    public boolean clear()
	{
		if(byteBuffer == null) return false;
		if(locked)
		{
			StsException.systemError("StsFileBlock.clear() shouldn't be called if locked is true. Block: " + debugString());
			return false;
		}
		try
		{
            if(StsFileBlocks.debug) StsException.systemDebug(this, "clear", " for file: " + fileMapBlocks.pathname + " block: " + nBlock);
            byteBuffer.force();
			byteBuffer.clear();
			StsToolkit.clean(byteBuffer);
			byteBuffer = null;
			if(fileMapBlocks.blockDeleted(size))return true;
			StsException.systemError("StsSeismicVolume.Plane.free() failed. " + nBlock + " already free.");
			return false;
		}
		catch(Exception e)
		{
			StsException.systemError(this, "clear", debugString());
            byteBuffer = null;
            return false;
		}
	}

	public long getMemoryUsed()
	{
		return size;
	}


	public String debugString()
	{
		if(fileMapBlocks == null)
			return new String("nullCube[" + nBlock + "]. locked: " + locked);
		else
			return new String(fileMapBlocks.pathname + "[" + nBlock + "]. locked: " + locked);
	}
}
