package com.Sts.Utilities.DataCube;

import com.Sts.MVC.*;
import com.Sts.Utilities.*;

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

/** StsBlocksMemoryManager is a singleton (currently in StsProject) and a static in StsFileBlocks which manages all
 *  StsFileBlocks volumes from StsSeismicVolume, StsVirtualVolume, StsPreStackLineSet3d, and StsVolumeConstructor
 *  (and any others).  When we wish to delete a block, we call this with the fileMappedBlocks from the calling
 *  class and the block in that blocks collection to be deleted.  These are matched up here, the block is deleted, and
 *  then the fileMappedBlocks itself is called to remove the block there.  Don't call StsFileMappedBlocks
 *  directly to do this; always call StsMappedBlocksMemory!
 */

public class StsBlocksMemoryManager
{
	long maxMemory;
	long memoryUsed = 0;
	double defaultMinFraction = 0.05;
	double defaultMaxFraction = 0.25;
	double minFraction = 0.1;
	double maxFraction = 0.5;
	int incSize = 1000;
	StsFileBlock[] blocks = new StsFileBlock[incSize]; // planes in memory
	public int nBlocksInMemory = 0;
	/** level 0: no mainDebug; level 1: summary mainDebug; level 2: detailed mainDebug */
	static final int debugLevel = 0;

	public StsBlocksMemoryManager(long maxMemoryMB)
	{
		this.maxMemory = maxMemoryMB * 1000000;
	}

	public void checkMemory()
	{
		if(memoryUsed < maxFraction * maxMemory)
		{
			return;
		}
		long memoryNeeded = memoryUsed - (long)(minFraction * maxMemory);
       
		long memoryFreed = 0;
		int nFreed = 0;
//        System.out.println("Number of planes in memory:" + nPlanesInMemory);
//        if(debugLevel > 1) System.out.println("nBlocks in memory: " + nBlocksInMemory);

        int[] deletableBlocks = new int[nBlocksInMemory];
		int nDeletableBlocks = 0;
		for(int n = 0; n < nBlocksInMemory; n++)
		{
			if(blocks[n] == null)
			{
				StsException.systemError(this, "checkMemory", " Block " + n + "  shouldn't be null");
				continue;
			}
			StsFileBlock block = blocks[n];
			if(block.byteBuffer != null && !block.locked)
			{
//				if(debugLevel > 1)StsException.systemDebug(this, "checkMemory", "deleting Block " + getDebugString(block, n) + " from mappedBlocksMemory ");
				deletableBlocks[nDeletableBlocks++] = n;
				long bytesFreed = block.size;
				memoryFreed += bytesFreed;
				if(memoryFreed > memoryNeeded)break;
			}
		}
		if(nDeletableBlocks == 0)return;

		removeDeletableBlocks(deletableBlocks, nDeletableBlocks);

		if(debugLevel > 0)
		{
			StsException.systemDebug(this, "checkMemory", "Memory freed. Blocks freed: " + nDeletableBlocks + " memoryFreed: " + memoryFreed + ". Current memoryUsed: " + memoryUsed);
        }
	}

	/** remove these blocks from dataCubeMemory. Calling routine is responsible for
	 *  actually deleting them from the dataPlanes array, however.
	 */
	/*
	 private int deleteBlocks(int nBlockMin, int nBlockMax)
	 {
	  StsFileBlock[] newBlocks = new StsFileBlock[blocks.length];
	  int nn = nBlockMin;
	  int nBlocksDeleted = 0;
	  int nBlocksAboveToMove = nBlocksInMemory-nBlockMax-1;
	  for(int n = nBlockMin; n <= nBlockMax; n++)
	  {
	   if(deleteBlock(n))
		nBlocksDeleted++;
	   else
		newBlocks[nn++] = blocks[n];
	  }
	  if(nBlocksDeleted == 0) return 0;
	  System.arraycopy(blocks, 0, newBlocks, 0, nBlockMin);
	  System.arraycopy(blocks, nBlockMax+1, newBlocks, nn, nBlocksAboveToMove);
	  blocks = newBlocks;
	  return nBlocksDeleted;
	 }
	 */
	public void printMemorySummary()
	{
		StsMessageFiles.infoMessage("Seismic planes in memory: " + nBlocksInMemory + ". Memory used: " + memoryUsed);
	}

	public void setMemoryFractions(float minFraction, float maxFraction)
	{
		this.minFraction = minFraction;
		this.maxFraction = maxFraction;
	}

	public void restoreMemoryFractions()
	{
		minFraction = defaultMinFraction;
		maxFraction = defaultMaxFraction;
	}

    public boolean clearBlock(StsFileBlock block)
    {
        if(block == null) return false;
        int index = StsMath.arrayGetIndex(blocks, block);
        if(index == StsParameters.NO_MATCH) return false;
        block.clear();
        StsMath.arrayDeleteElement(blocks, index);
        return true;
    }

    public boolean unlockClearBlock(StsFileBlock block)
    {
        if(block == null) return false;
        int index = StsMath.arrayGetIndex(blocks, block);
        if(index == StsParameters.NO_MATCH) return false;
        block.unlockClear();
        StsMath.arrayDeleteElement(blocks, index);
        return true;
    }
    
    public boolean clearBlocks(StsFileBlocks fileMapBlocks)
	{
		if(nBlocksInMemory == 0)return false;
		int[] deletableBlocks = new int[nBlocksInMemory];
		int nDeletableBlocks = 0;
		if(debugLevel > 0)System.out.println("removing all unlocked blocks in " + fileMapBlocks.pathname);

		for(int n = 0; n < nBlocksInMemory; n++)
		{
			if(blocks[n] == null)
			{
				StsException.systemError(this, "clearBlocks", " block shouldn't be null " + debugString(n));
				continue;
			}
			if(!blocks[n].locked && blocks[n].fileMapBlocks == fileMapBlocks)
			{
				deletableBlocks[nDeletableBlocks++] = n;
			}
		}
		if(nDeletableBlocks == 0)return true;
		return removeDeletableBlocks(deletableBlocks, nDeletableBlocks);
	}

    public boolean unlockAndClearBlocks(StsFileBlocks fileMapBlocks)
	{
		if(nBlocksInMemory == 0)return false;
		int[] deletableBlocks = new int[nBlocksInMemory];
		int nDeletableBlocks = 0;
		if(debugLevel > 0)System.out.println("removing all unlocked blocks in " + fileMapBlocks.pathname);

		for(int n = 0; n < nBlocksInMemory; n++)
		{
			if(blocks[n] == null)
			{
                StsException.systemError(this, "unlockAndClearBlocks", " Block shouldn't be null " + debugString(n));
				continue;
			}
			if(blocks[n].fileMapBlocks == fileMapBlocks)
			{
                blocks[n].locked = false;
                deletableBlocks[nDeletableBlocks++] = n;
			}
		}
		if(nDeletableBlocks == 0)return true;
		return removeDeletableBlocks(deletableBlocks, nDeletableBlocks);
	}

    public boolean clearUnlockedBlocks()
	{
		if(nBlocksInMemory == 0)return false;
		int[] deletableBlocks = new int[nBlocksInMemory];
		int nDeletableBlocks = 0;
		if(debugLevel > 0)System.out.println("removing all unlocked blocks");

		for(int n = 0; n < nBlocksInMemory; n++)
		{
			if(blocks[n] == null)
			{
				StsException.systemError(this, "clearUnlockedBlocks", " Block shouldn't be null " + debugString(n));
				continue;
			}
			if(!blocks[n].locked)
			{
				deletableBlocks[nDeletableBlocks++] = n;
			}
		}
		if(nDeletableBlocks == 0)return true;
		return removeDeletableBlocks(deletableBlocks, nDeletableBlocks);
	}

    private boolean removeDeletableBlocks(int[] deletableBlocks, int nDeletableBlocks)
	{
        if(nBlocksInMemory == 0) return true;
        StsFileBlock[] oldBlocks = new StsFileBlock[nBlocksInMemory];
        System.arraycopy(blocks, 0, oldBlocks, 0, nBlocksInMemory);
        int nextDeletableBlock = -1;
		int nNew = 0;
		int nOldBlocksInMemory = nBlocksInMemory;
		for(int n = 0; n < nDeletableBlocks; n++)
		{
			int prevDeletableBlock = nextDeletableBlock;
			nextDeletableBlock = deletableBlocks[n];
			for(int i = prevDeletableBlock + 1; i < nextDeletableBlock; i++)
			{
				blocks[nNew++] = oldBlocks[i];
			}
			clearBlock(nextDeletableBlock);
		}
		for(int i = nextDeletableBlock + 1; i < nOldBlocksInMemory; i++)
		{
			blocks[nNew++] = oldBlocks[i];
		}
        for(; nNew < nOldBlocksInMemory; nNew++)
            blocks[nNew] = null;
        if(debugLevel > 0) System.out.println("Blocks remaining in memory: " + nBlocksInMemory);
		return true;
	}


    private boolean clearBlock(int nBlock)
	{
		StsFileBlock block = blocks[nBlock];
        if(debugLevel > 0) System.out.println("    removing block: " + block.debugString());
        if(block.locked)
		{
			StsException.systemError("StsFileBlock.clear() shouldn't be called if locked is true.");
			return false;
		}
        boolean clearedOK = block.clear();
        nBlocksInMemory--;
        memoryUsed -= block.size;
        blocks[nBlock] = null;
        return clearedOK;
    }
/*
    private boolean removeDeletableBlocks(int[] deletableBlocks, int nDeletableBlocks)
	{
		StsFileBlock[] newBlocks = new StsFileBlock[blocks.length];
		int nextDeletableBlock = -1;
		int nNew = 0;
		int nCurrentBlocksInMemory = nBlocksInMemory;
		for(int n = 0; n < nDeletableBlocks; n++)
		{
			int prevDeletableBlock = nextDeletableBlock;
			nextDeletableBlock = deletableBlocks[n];
			for(int i = prevDeletableBlock + 1; i < nextDeletableBlock; i++)
			{
				newBlocks[nNew++] = blocks[i];
			}
			deleteBlock(nextDeletableBlock);
		}
		for(int i = nextDeletableBlock + 1; i < nCurrentBlocksInMemory; i++)
		{
			newBlocks[nNew++] = blocks[i];
		}
		blocks = newBlocks;
		if(debugLevel > 0) System.out.println("Blocks remaining in memory: " + nBlocksInMemory);
		return true;
	}
*/
	public boolean clearOtherBlocks(StsFileBlocks fileMapBlocks)
	{
		if(nBlocksInMemory == 0)return false;
		int[] deletableBlocks = new int[nBlocksInMemory];
		int nDeletableBlocks = 0;
		if(debugLevel > 0)System.out.println("removing all unlocked blocks other than those in " + fileMapBlocks.pathname);

		for(int n = 0; n < nBlocksInMemory; n++)
		{
			if(blocks[n] == null)
			{
				StsException.systemError(this, "clearOtherBlocks(StsFileBlocks)", " Block shouldn't be null " + debugString(n));
				continue;
			}
			if(!blocks[n].locked && blocks[n].fileMapBlocks != fileMapBlocks)
			{
				deletableBlocks[nDeletableBlocks++] = n;
			}
		}
		if(nDeletableBlocks == 0)return true;
		return removeDeletableBlocks(deletableBlocks, nDeletableBlocks);
	}

	public boolean clearOtherBlocks(StsFileBlocks[] filesMapBlocks)
	{
		if(nBlocksInMemory == 0)return false;
		int[] deletableBlocks = new int[nBlocksInMemory];
		int nDeletableBlocks = 0;

		if(debugLevel > 0)
		{
			System.out.print("removing all unlocked blocks other than those in ");
			for(int i = 0; i < filesMapBlocks.length; i++)
				System.out.print(" " + filesMapBlocks[i].pathname);
			System.out.println("");
		}

		for(int n = 0; n < nBlocksInMemory; n++)
		{
			if(blocks[n] == null)
			{
				StsException.systemError(this, "clearOtherBlocks(StsFileBlocks[]", " Block shouldn't be null " + debugString(n));
				continue;
			}
			if(!blocks[n].locked && !blockIsInSets(blocks[n], filesMapBlocks))
				deletableBlocks[nDeletableBlocks++] = n;
		}
		if(nDeletableBlocks == 0)return true;
		return removeDeletableBlocks(deletableBlocks, nDeletableBlocks);
	}

	private boolean blockIsInSets(StsFileBlock block, StsFileBlocks[] filesMapBlocks)
	{
		for(int i = 0; i < filesMapBlocks.length; i++)
			if(block.fileMapBlocks == filesMapBlocks[i])return true;
		return false;
	}

	public boolean clearBlocks(StsFileBlocks[] filesMapBlocks)
	{
		if(nBlocksInMemory == 0)return false;
		int[] deletableBlocks = new int[nBlocksInMemory];
		int nDeletableBlocks = 0;
		for(int n = 0; n < nBlocksInMemory; n++)
		{
			if(blocks[n] == null)
			{
				StsException.systemError(this, "clearBlocks(StsFileBlocks[])", " Block shouldn't be null " + debugString(n));
				continue;
			}
			for(int i = 0; i < filesMapBlocks.length; i++)
			{
				StsFileBlocks fileMapBlocks = filesMapBlocks[i];
				if(debugLevel > 0) System.out.println("removing all unlocked blocks in " + fileMapBlocks.pathname);
				if(blocks[n] == null)
				{
					StsException.systemError("StsMappedBlocksMemory.clearBlocks(filesMapBlocks) failed for block " + n);
					continue;
				}
				if(!blocks[n].locked && blocks[n].fileMapBlocks != fileMapBlocks)
				{
					deletableBlocks[nDeletableBlocks++] = n;;
					break;
				}
			}
		}
		if(nDeletableBlocks == 0) return true;
		return removeDeletableBlocks(deletableBlocks, nDeletableBlocks);
	}

    /** removes all blocks from memory.  Calling routine is responsible for removing these same blocks from any block
     *  collection (subclasses of StsFileBlocks).
     * @return
     */
    public boolean clearAllBlocks()
	{
		if(nBlocksInMemory == 0) return false;
		if(debugLevel > 0) System.out.println("removing all " + nBlocksInMemory + " blocks in memory ");
		for(int n = 0; n < nBlocksInMemory; n++)
		{
            if(blocks[n] == null)
            {
                StsException.systemError(this, "clearAllBlocks", " Block shouldn't be null " + debugString(n));
                continue;
            }
            blocks[n].locked = false;
			blocks[n].clear();
            blocks[n] = null;
        }
		memoryUsed = 0;
		nBlocksInMemory = 0;
		return true;
	}

    public void addBlock(StsFileBlock block)
	{
        if(block == null)
        {
            StsException.systemError(this, "addBlock", "Attempting to add a null block.");
            return;
        }
        blocks = (StsFileBlock[])StsMath.arrayAddElement(blocks, block, nBlocksInMemory, incSize);
        int location = nBlocksInMemory;
        nBlocksInMemory++;
		this.memoryUsed += block.size;
		if(debugLevel > 1)System.out.println("added Block " + getDebugString(block, location));
	}

	private String getDebugString(StsFileBlock block, int location)
	{
		if(block == null) return "null mappedBlockObject";
		else return block.debugString() + " " + debugString(location);
	}

	public String debugString(int location)
	{
		return " blocksMemoryManager location: " + location + " Blocks in memory:" + nBlocksInMemory + " MemoryUsed: " + memoryUsed;
	}


	/*
	 public boolean removePlane(StsDataPlane dataPlane)
	 {
	  if (dataPlane == null || dataPlanes == null) return false;
	  StsDataCube dataCube = dataPlane.dataCube;
	  int nPlane = dataPlane.nPlane;
	  int planeMemory = dataCube.getPlaneMemorySize();

	  for (int n = 0; n < nPlanesInMemory; n++)
	  {
	   if (dataPlanes[n] == null)
	   {
		StsException.systemError("StsSeismicVolume.dataCube.removePlane() failed for " + getDebugString(dataCube, nPlane));
		continue;
	   }
	   if (dataPlanes[n].dataCube == dataCube && dataPlanes[n].nPlane == nPlane)
	   {
		dataPlanes[n].free();
		dataPlanes = (StsDataPlane[]) StsMath.arrayDeleteElement(dataPlanes, n);
		nPlanesInMemory--;
		memoryUsed -= planeMemory;
		if(mainDebug) System.out.println("removed Plane from dataCubeMemory " + getDebugString(dataCube, nPlane));
		return true;
	   }
	  }
	  return false;
	 }
	 */

	// remove blocks not between keepMin and keepMax as long as they are not between excludeMin and excludeMax
/*
	public void clearBlocks(StsCubeFileBlocks fileMapBlocks, int keepMin, int keepMax, int excludeMin, int excludeMax)
	{
		if(blocks == null)return;

		for(int n = nBlocksInMemory - 1; n >= 0; n--)
		{
			if(blocks[n] == null)
			{
				StsException.systemError("StsMappedBlockMemory.removePlanes() failed for " + getDebugString(fileMapBlocks, -1));
				continue;
			}
			if(blocks[n].fileMapBlocks == fileMapBlocks)
			{
				int nBlock = blocks[n].nBlock;
				if((nBlock < keepMin || nBlock > keepMax) && (nBlock < excludeMin || nBlock > excludeMax))
				{
					deleteBlock(n);
					blocks = (StsFileBlock[])StsMath.arrayDeleteElement(blocks, n);
				}
			}
		}
	}
*/
	/*
	 public void removeAll()
	 {
	  if (dataPlanes == null)
	  {
	   return;
	  }
	  for (int n = 0; n < nPlanesInMemory; n++)
	  {
	   dataPlanes[n] = null;
	  }
	  nPlanesInMemory = 0;
	  memoryUsed = 0;
	 }
	 */
}
