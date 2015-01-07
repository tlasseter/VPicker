package com.Sts.Utilities.DataCube;

import com.Sts.DBTypes.*;
import com.Sts.Utilities.*;

import java.io.*;

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
public class StsPreStackLineFileBlocks extends StsFileBlocks
{
	int nGathers;
	int nSlices;
    static StsTimer openFileTimer;
    static StsTimer readFileTimer;
	static final boolean debug = false;

    /** scratch float data filled in accessing all traces in gather */
    static float[] floatData = null;

    static
	{
	    if(debug)
		{
			openFileTimer = StsTimer.constructor("open/close prestack files timer");
			readFileTimer = StsTimer.constructor("read prestack files timer");
		}
	}
    public StsPreStackLineFileBlocks(int nGathers, int maxGatherTraces, int[] nLastTraceInGathers, int nSlices, StsPreStackLineSet lineSet, String filename, int nBytesPerSample,
									 StsBlocksMemoryManager blocksMemoryManager) throws FileNotFoundException, IOException
	{
		this.nGathers = nGathers;
		this.nSlices = nSlices;
        int nTotalTraces = nLastTraceInGathers[nGathers - 1] + 1;
		fileSize = nTotalTraces * nSlices * nBytesPerSample;

        initialize(lineSet, filename, "r", nBytesPerSample, blocksMemoryManager);
//		approxBlockSize = 1000000;
//		maxBlockSize = 2000000;
// Moved to getGatherDataf to solve too many files open problem -----	classInitialize(directory, filename, "rw", nBytesPerSample, blocksMemoryManager);

		long maxBytesPerGather = maxGatherTraces * nSlices * nBytesPerSample; // approx size of a gather in bytes

        if(fileSize < maxBlockSize)
		{
			nTotalBlocks = 1;
			blocks = new StsFileBlock[1];
			blocks[0] = new StsPreStackLineFileBlock(this, 0, 0, fileSize, 0, nGathers - 1, 0, nTotalTraces - 1);
		}
		else
		{
            int nGathersPerBlock = (int)(maxBlockSize/maxBytesPerGather);
//			int nGathersPerBlock = Math.round(approxBlockSize / maxBytesPerGather);
			nGathersPerBlock = Math.max(1, nGathersPerBlock);
			nTotalBlocks = StsMath.ceiling((float)nGathers / nGathersPerBlock);
			long nBytesPerBlock = fileSize / nTotalBlocks;
			blocks = new StsFileBlock[nTotalBlocks];
			long offset = 0;
			long nBlockBytes = 0;
			int nLastTrace = -1;
			int nLastBlockTrace = -1;
			int nLastBlockGather = -1;
			int nBlock = 0;
			for(int n = 0; n < nGathers; n++)
			{
				int nFirstTrace = nLastTrace + 1;
				nLastTrace = nLastTraceInGathers[n];
				int nTraces = nLastTrace - nFirstTrace + 1;
				int nGatherBytes = nTraces * nSlices * nBytesPerSample;
				nBlockBytes += nGatherBytes;
				if(nBlockBytes >= nBytesPerBlock)
				{
					int nFirstBlockTrace = nLastBlockTrace + 1;
					int nFirstBlockGather = nLastBlockGather + 1;
					nLastBlockTrace = nLastTrace;
					nLastBlockGather = n;
					blocks[nBlock] = new StsPreStackLineFileBlock(this, nBlock, offset, nBlockBytes, nFirstBlockGather, nLastBlockGather, nFirstBlockTrace, nLastBlockTrace);
					nBlock++;
					offset += nBlockBytes;
					nBlockBytes = 0;
				}
			}
			if(nBlockBytes > 0)
			{
				int nFirstBlockTrace = nLastBlockTrace + 1;
				nLastBlockTrace = nLastTrace;
				int nFirstBlockGather = nLastBlockGather + 1;
				nLastBlockGather = nGathers - 1;
				blocks[nBlock] = new StsPreStackLineFileBlock(this, nBlock, offset, nBlockBytes, nFirstBlockGather, nLastBlockGather, nFirstBlockTrace, nLastBlockTrace);
				nBlock++;
			}
			nTotalBlocks = nBlock;
			if(nTotalBlocks < blocks.length)
				blocks = (StsFileBlock[])StsMath.trimArray(blocks, nTotalBlocks);
		}
	}
/*
	public FloatBuffer getGatherData(int nGather, int nFirstTrace, int nLastTrace, int nCroppedSlices, boolean applyAGC, int windowWidth, boolean AGCchanged)
	{
		if(nGather < 0 || nGather >= nGathers)return null;
		for(int n = 0; n < nTotalBlocks; n++)
		{
			StsPreStackLineFileBlock block = (StsPreStackLineFileBlock)blocks[n];
			if(nFirstTrace >= block.nFirstBlockTrace && nLastTrace <= block.nLastBlockTrace)
			{
				return block.getGatherData(nGather, nFirstTrace, nLastTrace, nCroppedSlices, applyAGC, windowWidth, AGCchanged);
			}
		}
		return null;
	}
*/
	public float[] getGatherDataf(StsBlocksMemoryManager blocksMemoryManager, int nBytesPerSample, int nGather, int nFirstTrace, int nLastTrace, int nSlices)
	{
        // Moved from constructor because we have too many files open, need to open/read/close everytime.
        if(debug) openFileTimer.start();
//        initialize(pathname, "rw", nBytesPerSample, blocksMemoryManager);
        if(debug) openFileTimer.stopAccumulateIncrementCount();
		if(nGather < 0 || nGather >= nGathers)return null;
		for(int n = 0; n < nTotalBlocks; n++)
		{
			StsPreStackLineFileBlock block = (StsPreStackLineFileBlock)blocks[n];
			if(nFirstTrace >= block.nFirstBlockTrace && nLastTrace <= block.nLastBlockTrace)
			{
                if(debug) readFileTimer.start();
                float[] data = block.getGatherDataf(nGather, nFirstTrace, nLastTrace, nSlices);
                if(debug)
				{
					readFileTimer.stopAccumulateIncrementCount();
					openFileTimer.start();
				}
//                closeFile();
				if(debug)
				{
					openFileTimer.stopAccumulateIncrementCount();
					openFileTimer.printElapsedTime();
					readFileTimer.printElapsedTime();
				}
                return data;
			}
		}
		return null;
	}

    static public float[] getScratchFloatData(int nFloats)
    {
        if(floatData == null || floatData.length < nFloats)
            floatData = new float[nFloats];
        return floatData;
    }
}
