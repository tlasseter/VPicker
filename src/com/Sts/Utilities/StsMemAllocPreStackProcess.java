package com.Sts.Utilities;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsMemAllocPreStackProcess extends StsMemoryAllocation
{
	int nTotalTraces;
    int nBytesPerOutputTrace;
	int maxInputMegaBytes = 20; // amount of page memory used by a single file NIO mapped byte buffer which is the input SEGY file
	int nBytesPerOutputSample = 1;
	int bytesPerInputTrace;
	int nSamplesPerTrace;
	public int nInputBytesPerBlock;
//	public int nOutputSamplesPerWrite;
	public int nTracesPerBlock;
	public int nTotalBytes;

//	public byte[] blockOffsetBytes;
	public int nBlocks;
	public int nBlocksPerWrite;

	static public boolean dontAskAgain = false;

	private StsMemAllocPreStackProcess(int nTotalTraces, int bytesPerInputTrace, int nSamplesPerTrace) throws StsException
	{
        initializeMemoryAllocation();
        this.nTotalTraces = nTotalTraces;
		this.bytesPerInputTrace = bytesPerInputTrace;
        this.nSamplesPerTrace = nSamplesPerTrace;
        nTracesPerBlock = (int) (maxInputMegaBytes * 1000000 / bytesPerInputTrace);
		nBlocks = StsMath.ceiling(((float) nTotalTraces) / nTracesPerBlock);
		nInputBytesPerBlock = nTracesPerBlock * bytesPerInputTrace;
		nBytesPerOutputTrace = nSamplesPerTrace * nBytesPerOutputSample; // we are writing offset-ordered output
		nBlocksPerWrite = 1;
    }

    static public StsMemAllocPreStackProcess constructor(int nTotalTraces, int bytesPerInputTrace, int nSamplesPerTrace)
	{
		try
		{
			return new StsMemAllocPreStackProcess(nTotalTraces, bytesPerInputTrace, nSamplesPerTrace);
		}
		catch(Exception e)
		{
			StsException.systemError(e.getMessage());
			return null;
		}
	}

    protected boolean nextIterationOK() { return true; }

    protected boolean allocateMemoryIteration()
	{
        return true;
	}

	public void setDontAskAgain(boolean b)
	{
		dontAskAgain = b;
	}

	public boolean getDontAskAgain()
	{
		return dontAskAgain;
	}

    public void freeMemory() { }
}
