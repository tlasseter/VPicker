package com.Sts.Utilities.DataCube;

import com.Sts.DBTypes.*;
import com.Sts.Utilities.*;

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
public class StsPreStackLineFileBlock extends StsFileBlock
{
	int nFirstGather;
	int nLastGather;
	int nFirstBlockTrace;
	int nLastBlockTrace;
	byte[] gatherStatus;

    static final byte NOT_IN_MEMORY = 0;
	static final byte IN_MEMORY_AND_SCALED = 1;
	static final byte AGC_APPLIED = 2;

	public StsPreStackLineFileBlock(StsFileBlocks fileMapBlocks, int nBlock, long offset, long size, int nFirstGather, int nLastGather, int nFirstBlockTrace, int nLastBlockTrace)
	{
		super(fileMapBlocks, nBlock, offset, size);
		this.nFirstBlockTrace = nFirstBlockTrace;
		this.nLastBlockTrace = nLastBlockTrace;
		this.nFirstGather = nFirstGather;
		this.nLastGather = nLastGather;
		int nBlockGathers = nLastGather - nFirstGather + 1;
		gatherStatus = new byte[nBlockGathers];
	}

/*
	  public FloatBuffer getGatherData(int nGather, int nFirstTrace, int nLastTrace, int nCroppedSlices, boolean applyAGC, int windowWidth, boolean AGCchanged)
	{
		FloatBuffer floatBuffer;
		int blockOffset = (nFirstTrace - nFirstBlockTrace)*nCroppedSlices*4;
		byte status = getStatus(nGather);

		if ((status == AGC_APPLIED && !applyAGC) || AGCchanged)
		{
			clear();
		}

		if (status == NOT_IN_MEMORY || (status == AGC_APPLIED && !applyAGC) || (status != AGC_APPLIED && applyAGC) || AGCchanged)
		{
			int nGatherTraces = nLastTrace - nFirstTrace + 1;
			int nFloats = nGatherTraces*nCroppedSlices;
			float [] data = new float[nFloats];
			if(StsPreStackLineSet3d.debugTimer) StsSeismicTimer.getInputBlockTimer.start();
			MappedByteBuffer byteBuffer = fileMapBlocks.getBlock(nBlock, FileChannel.MapMode.READ_ONLY);
			if(StsPreStackLineSet3d.debugTimer) StsSeismicTimer.getInputBlockTimer.stop();
			if(StsPreStackLineSet3d.debugTimer) StsSeismicTimer.readInputBlockTimer.start();
			{
				byteBuffer.position(blockOffset);
			    floatBuffer = byteBuffer.asFloatBuffer().slice();
				floatBuffer.get(data);
				floatBuffer.rewind();
			}
			if(StsPreStackLineSet3d.debugTimer) StsSeismicTimer.readInputBlockTimer.stop();

			if(StsPreStackLineSet3d.debugTimer) StsSeismicTimer.normalizeTimer.start();
			StsSeismicFilter.normalizeAmplitude(floatBuffer, data, nFirstTrace, nLastTrace, nCroppedSlices);
			if(StsPreStackLineSet3d.debugTimer) StsSeismicTimer.normalizeTimer.stop();
			setStatus(nGather, IN_MEMORY_AND_SCALED);
			status = IN_MEMORY_AND_SCALED;
		}
		else
		{
			byteBuffer.position(blockOffset);
			floatBuffer = byteBuffer.asFloatBuffer().slice();
		}

		if(status == IN_MEMORY_AND_SCALED && applyAGC)
		{
			StsSeismicFilter.applyAGC(floatBuffer, nFirstTrace, nLastTrace, nCroppedSlices, windowWidth);
			int nGatherTraces = nLastTrace - nFirstTrace + 1;
            int nFloats = nGatherTraces*nCroppedSlices;
            float [] data = new float[nFloats];

			floatBuffer.get(data);
			floatBuffer.rewind();

			StsSeismicFilter.normalizeAmplitude(floatBuffer, data, nFirstTrace, nLastTrace, nCroppedSlices);
			setStatus(nGather, AGC_APPLIED);
		}
		return floatBuffer;
	}
	*/

	/*
	public float[] getGatherDataf(int nGather, int nFirstTrace, int nLastTrace, int nCroppedSlices, boolean applyAGC, int windowWidth, boolean AGCchanged)
	{
		FloatBuffer floatBuffer;
		int blockOffset = (nFirstTrace - nFirstBlockTrace)*nCroppedSlices*4;
		byte status = getStatus(nGather);
		int nGatherTraces = nLastTrace - nFirstTrace + 1;
		int nFloats = nGatherTraces*nCroppedSlices;
		float [] data=null;

// in this formulation we don't care about "status" -- we always fetch up the data from a buffer
		if (status == NOT_IN_MEMORY )
		{
			if(StsPreStackLineSet3d.debugTimer) StsSeismicTimer.getInputBlockTimer.start();
			MappedByteBuffer byteBuffer = fileMapBlocks.getBlock(nBlock, FileChannel.MapMode.READ_ONLY);
			if(StsPreStackLineSet3d.debugTimer) StsSeismicTimer.getInputBlockTimer.stop();

	        if(StsPreStackLineSet3d.debugTimer) StsSeismicTimer.readInputBlockTimer.start();
			{
				byteBuffer.position(blockOffset);
				floatBuffer = byteBuffer.asFloatBuffer().slice();
				data = new float[nFloats];
				floatBuffer.get(data);
				floatBuffer.rewind();
			}
			if(StsPreStackLineSet3d.debugTimer) StsSeismicTimer.readInputBlockTimer.stop();

			setStatus(nGather, IN_MEMORY_AND_SCALED);
            status = IN_MEMORY_AND_SCALED;
		}
		else
		{
			byteBuffer.position(blockOffset);
			floatBuffer = byteBuffer.asFloatBuffer().slice();
			data = new float[nFloats];
			floatBuffer.get(data);
			floatBuffer.rewind();
		}

		// alway have to nromalize & agc as the floatbuffer is read-only here

		if(StsPreStackLineSet3d.debugTimer) StsSeismicTimer.normalizeTimer.start();
	    StsSeismicFilter.normalizeAmplitude(data, nFirstTrace, nLastTrace, nCroppedSlices);
	    if(StsPreStackLineSet3d.debugTimer) StsSeismicTimer.normalizeTimer.stop();

		if(applyAGC)
		{
			StsSeismicFilter.applyAGC(data, nFirstTrace, nLastTrace, nCroppedSlices, windowWidth);
			StsSeismicFilter.normalizeAmplitude(data, nFirstTrace, nLastTrace, nCroppedSlices);
			setStatus(nGather, AGC_APPLIED);
		}


		return data;
}
	*/
	public float[] getGatherDataf(int nGather, int nFirstTrace, int nLastTrace, int nSlices)
	{
		int blockOffset = (nFirstTrace - nFirstBlockTrace)*nSlices*4;
		byte status = getStatus(nGather);
		int nGatherTraces = nLastTrace - nFirstTrace + 1;
		int nFloats = nGatherTraces*nSlices;

// in this formulation we don't care about "status" -- we always fetch up the data from a buffer
		if (status == NOT_IN_MEMORY )
		{
			if(StsPreStackLineSet3d.debugTimer) StsSeismicTimer.getInputBlockTimer.start();
			MappedByteBuffer byteBuffer = fileMapBlocks.getBlockBuffer(nBlock, FileChannel.MapMode.READ_ONLY);
            if(byteBuffer == null) return null;
            if(StsPreStackLineSet3d.debugTimer) StsSeismicTimer.getInputBlockTimer.stopAccumulateIncrementCount();
            setStatus(nGather, IN_MEMORY_AND_SCALED);
			status = IN_MEMORY_AND_SCALED;
		}
        return readFloatData(byteBuffer, blockOffset, nFloats);
    }

    private float[] readFloatData(MappedByteBuffer byteBuffer, int blockOffset, int nFloats)
    {
        if(StsPreStackLineSet3d.debugTimer) StsSeismicTimer.readInputBlockTimer.start();

        byteBuffer.position(blockOffset);
        FloatBuffer floatBuffer = byteBuffer.asFloatBuffer();
        float[] floatData = StsPreStackLineFileBlocks.getScratchFloatData(nFloats);
        floatBuffer.get(floatData, 0, nFloats);
        floatBuffer.rewind();

        if(StsPreStackLineSet3d.debugTimer) StsSeismicTimer.readInputBlockTimer.stopAccumulateIncrementCount();
        return floatData;
    }

    public boolean clear()
	{
		for(int n = 0; n < gatherStatus.length; n++)
		{
			gatherStatus[n] = NOT_IN_MEMORY;
		}
		return super.clear();
	}

	private byte getStatus(int nGather)
	{
		return gatherStatus[nGather - nFirstGather];
	}

	private void setStatus(int nGather, byte status)
	{
		gatherStatus[nGather - nFirstGather] = status;
	}


}
