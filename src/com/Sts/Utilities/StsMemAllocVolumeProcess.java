package com.Sts.Utilities;

import com.Sts.DBTypes.*;
import com.Sts.Types.*;

/** This class allocates 3 chunks of memory used for byte volume output: outputRowPlaneBytes, outputColBlockBytes, outputSliceBlockBytes
 *  The size of the first is fixed (nInputCols*nSamples), but the other two are sized by nRowsPerBlock*nInputCols*nSamples.
 *  Allocated a specified fraction of free memory, we allocate as many rowsPerBlock as can be accommodated.
 */
public class StsMemAllocVolumeProcess extends StsMemoryAllocation
{
    public int nInputRows;
    public int nInputCols;
    public int nInputSlices;
    public int nOutputRows;
    public int nOutputCols;
    public int nOutputSlices;
    public int nBytesPerInputSample;
    public int nBytesPerOutputSample;
    public int nInputVolumes = 1;
    public long nInputBytesPerRow;
//    public int nOutputBytesPerRow;
    public int nInputRowsPerBlock;
    public int nOutputRowsPerBlock;
    public int nInputBlocks;
    public int nOutputBlocks;
    public int nInputBytesPerBlock;
//    public int nBlocksPerWrite;
    public int nOutputBytesPerBlock;
    public int nInputBlockTraces;
    public int nOutputSamplesPerRow;
    public int nInputSamplesPerBlock;
    public int nOutputSamplesPerInputBlock;
    public byte[] outputRowPlaneBytes;
	public byte[] outputColBlockBytes;
	public byte[] outputSliceBlockBytes;

    private StsMemAllocVolumeProcess(StsSeismicVolume inputVolume, StsSeismicVolume outputVolume) throws StsException
    {
        initializeMemoryAllocation();
        setupInputMemory(inputVolume);
        setupOutputMemory(outputVolume, false);
    }

    private StsMemAllocVolumeProcess(StsSegyVolume inputVolume, StsCroppedBoundingBox cropBox, boolean useOutputMappedBuffers) throws StsException
    {
        initializeMemoryAllocation();
        setupInputMemory(inputVolume);
        setupOutputMemory(cropBox, useOutputMappedBuffers);
    }

    private StsMemAllocVolumeProcess(StsSeismicVolume outputVolume) throws StsException
    {
        initializeMemoryAllocation();
        setupInputMemory(outputVolume);
        setupOutputMemory(outputVolume, false);
    }

    private void setupInputMemory(StsSeismicVolume volume)
    {
        nInputRows = volume.nRows;
        nInputCols = volume.nCols;
        nInputSlices = volume.nSlices;
        if(volume.isDataFloat)
        {
            nBytesPerInputSample = 4;
            nInputBlocks = volume.fileMapRowFloatBlocks.nTotalBlocks;
            nInputRowsPerBlock = (int)volume.fileMapRowFloatBlocks.nPlanesPerBlock;
        }
        else
        {
            nBytesPerInputSample = 1;
            nInputBlocks = volume.filesMapBlocks[1].nTotalBlocks;
            nInputRowsPerBlock = (int)volume.filesMapBlocks[1].nPlanesPerBlock;
        }
        nInputBytesPerRow = nInputCols*nInputSlices*nBytesPerInputSample;
        nInputBlockTraces = nInputRowsPerBlock*nInputCols;
        nInputSamplesPerBlock = nInputBlockTraces*nInputSlices;
        nInputBytesPerBlock = nInputSamplesPerBlock*nBytesPerInputSample;
    }

    private void setupInputMemory(StsSegyVolume segyVolume)
    {
        // we need to read all the rows from the first to last row of the cropped volume even though some may be decimated
        if(segyVolume.cropBox.isCropped)
            nInputRows = segyVolume.cropBox.rowMax - segyVolume.cropBox.rowMin + 1;
        else
            nInputRows = segyVolume.nRows;

        nInputCols = segyVolume.nCols;
        nInputSlices = segyVolume.nSlices;

        nBytesPerInputSample = segyVolume.segyData.bytesPerSample;
        int nBytesPerTrace = segyVolume.segyData.bytesPerTrace;
        // If volume is regular, we can compute an exact number of bytesPerRow and bytes per block.
        // If irregular, we can only approximate the number of bytes per block from an estimated rows per block.
        // For the irregular case, we have to search for last trace header in block and end of last block row,
        // compute the number of traces per block and bytes per block allowing us to map the byteBuffer to the block.
        nInputBytesPerRow = nInputCols*nBytesPerTrace;
        nInputRowsPerBlock = (int)Math.max(1, maxMapBufferSize/nInputBytesPerRow);
        nInputRowsPerBlock = Math.min(nInputRows, nInputRowsPerBlock);
        nInputBlocks = nInputRows/nInputRowsPerBlock;
        if(nInputBlocks*nInputRowsPerBlock < nInputRows) nInputBlocks++;
        nInputBlockTraces = nInputRowsPerBlock*nInputCols;
        if(nInputRowsPerBlock*nInputBlocks < nInputRows) nInputBlocks++;

        nInputBytesPerBlock = nInputBlockTraces*nBytesPerTrace;
        nInputSamplesPerBlock = nInputBlockTraces*nInputSlices;
    }

    private void setupOutputMemory(StsRotatedGridBoundingBox volume, boolean useOutputMappedBuffers) throws StsException
    {
        nOutputRows = volume.nRows;
        nOutputCols = volume.nCols;
        nOutputSlices = volume.nSlices;
        nBytesPerOutputSample = 4;

        nOutputSamplesPerRow = nOutputCols*nOutputSlices;
        // mapBuffer blockSize is controlled by the number of bytes in a float block we are processing
        // we will write out a rpo
        int nOutputBytesPerRow = nOutputSamplesPerRow*nBytesPerOutputSample;
        nOutputRowsPerBlock = Math.max(1, (int)(maxMemoryToUse-nOutputSamplesPerRow)/(2*nOutputSamplesPerRow));
    /*
        if(useOutputMappedBuffers)
            nOutputRowsPerBlock = Math.max(1, (int)(maxMapBufferSize/nOutputSamplesPerRow));
        else
            nOutputRowsPerBlock = Math.max(1, (int)(maxMemoryToUse-nOutputSamplesPerRow)/(2*nOutputSamplesPerRow));
    */
        nOutputRowsPerBlock = Math.min(nOutputRows, nOutputRowsPerBlock);
        // not currently used
        nOutputBytesPerBlock = nOutputRowsPerBlock*nOutputSamplesPerRow;
        nOutputBlocks = nOutputRows/nOutputRowsPerBlock;
        if(nOutputRowsPerBlock*nOutputBlocks < nOutputRows) nOutputBlocks++;
        nOutputSamplesPerInputBlock = nInputRowsPerBlock*nOutputSamplesPerRow;
    }

    static public StsMemAllocVolumeProcess constructor(StsSeismicVolume inputVolume, StsSeismicVolume outputVolume)
    {
        try
        {
            return new StsMemAllocVolumeProcess(inputVolume, outputVolume);
        }
 		catch(Exception e)
		{
			StsException.systemError(e.getMessage());
			return null;
		}
    }

    static public StsMemAllocVolumeProcess constructor(StsSegyVolume segyVolume, StsCroppedBoundingBox cropBox, boolean usingOutputMappedBuffers)
    {
        try
        {
            return new StsMemAllocVolumeProcess(segyVolume, cropBox, usingOutputMappedBuffers);
        }
 		catch(Exception e)
		{
			StsException.systemError(e.getMessage());
			return null;
		}
    }

    static public StsMemAllocVolumeProcess constructor(StsSeismicVolume outputVolume)
    {
        try
        {
            return new StsMemAllocVolumeProcess(outputVolume);
        }
 		catch(Exception e)
		{
			StsException.systemError(e.getMessage());
			return null;
		}
    }

    protected boolean nextIterationOK() { return nOutputRowsPerBlock > 1; }

    protected boolean allocateMemoryIteration()
	{
		int nTotalBytesPerOutputBlock = 0;
        try
		{
            outputRowPlaneBytes = new byte[nOutputSamplesPerRow];
            nOutputRowsPerBlock = (int)(maxMemoryToUse - nOutputSamplesPerRow)/(2* nOutputSamplesPerRow);
            nOutputRowsPerBlock = StsMath.minMax(nOutputRowsPerBlock, 1, nOutputRows);
            nOutputBytesPerBlock = nOutputRowsPerBlock * nOutputCols * nOutputSlices;
			outputColBlockBytes = new byte[nOutputBytesPerBlock];
			outputSliceBlockBytes = new byte[nOutputBytesPerBlock];
            nTotalBytesPerOutputBlock = nOutputSamplesPerRow + 2* nOutputBytesPerBlock;
            checkMemoryStatus("allocate memory succeeded for " + nOutputRowsPerBlock + " rowsPerBlock and " + nTotalBytesPerOutputBlock + " bytes");
		}
		catch (OutOfMemoryError e)
		{
			freeMemory();
			checkMemoryStatus("allocate memory failed for " + nOutputRowsPerBlock + " rowsPerBlock and " + nTotalBytesPerOutputBlock + " bytes");
			nOutputRowsPerBlock /= 2;
			return false;
		}
		return true;
	}

    public void freeMemory()
    {
        outputRowPlaneBytes = null;
        outputColBlockBytes = null;
	    outputColBlockBytes = null;
    }
}
