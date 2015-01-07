package com.Sts.SeismicAttributes;

import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

/**
 *
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002-2003</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version beta 1.0
 */

/** Given a set of inputVolumes required for constructing this seismic attribute, construct the attribute volume.
 *  Currently implementd for Hilbert attributes which are subclasses of this abstract class.
 */
public  class StsGeneralVolumeConstructor extends StsVolumeConstructor
{
    public StsSeismicBoundingBox seismicBoundingBox;
    boolean inputIsFloat = true;
	StsTimer timer = null;
	boolean runTimer = false;
	private int nSlices;
	String stsDirectory;
	String rowFloatFilename;

    static final float nullValue = StsParameters.nullValue;

	//static methods
    static public StsSeismicVolume createInterpolatedVolume(StsModel model, StsSeismicBoundingBox boundingBox, String volumeName, StsProgressPanel panel)
	{
		if(panel != null) panel.appendLine("Creating " + volumeName + " volume.");

		StsGeneralVolumeConstructor volumeConstructor = StsGeneralVolumeConstructor.constructor(model, boundingBox, panel);

		if(panel != null) panel.appendLine("Successfully created " + volumeName + " volume.");

		return volumeConstructor.outputVolume;
	}

    static public StsGeneralVolumeConstructor constructor(StsModel model, StsSeismicBoundingBox boundingBox, StsProgressPanel panel)
    {
        try
        {
            return new StsGeneralVolumeConstructor(model, boundingBox, panel);
        }
        catch (Exception e)
        {
            StsMessage.printMessage("StsGeneralVolumeConstructor.constructor() failed.");
            return null;
        }
    }

    //private constructor
    private StsGeneralVolumeConstructor(StsModel model, StsSeismicBoundingBox seismicBoundingBox, StsProgressPanel panel)
    {
        this.model = model;
        this.seismicBoundingBox = seismicBoundingBox;
        this.panel = panel;
        outputVolume = StsSeismicVolume.initializeAttributeVolume(model, seismicBoundingBox, dataMin, dataMax, true, true, seismicBoundingBox.getName(), "rw");
        createOutputVolume();
    }


    public StsGeneralVolumeConstructor()
	{
	}

    public void createVolume(StsSeismicBoundingBox inputVolume)
	{
		this.seismicBoundingBox = inputVolume;
        createOutputVolume();
    }

    public void initializeVolumeInput()
    {
        if(seismicBoundingBox == null) return;
        super.initializeVolumeInput();
        outputFloatBuffer = outputVolume.createMappedFloatRowBuffer();
        outputPosition = 0;
        traceFloats = null;
    }

    public boolean allocateMemory()
    {
        memoryAllocation = StsMemAllocVolumeProcess.constructor(outputVolume);
        if(memoryAllocation == null) return false;
        nOutputRowsPerBlock = memoryAllocation.nOutputRowsPerBlock;
        nOutputBlocks = memoryAllocation.nOutputBlocks;
        nOutputSamplesPerInputBlock = memoryAllocation.nOutputSamplesPerInputBlock;
        nOutputSamplesPerRow = memoryAllocation.nOutputSamplesPerRow;
        return true;
    }

	public  boolean processBlockInput(int nBlock, String mode ) { return processBlockInput(nBlock); }
    public  boolean processBlockInput(int nBlock)
    {
        return doProcessInputBlock(nBlock);
    }

    public boolean doProcessInputBlock(int nBlocks)
    {
        if (runTimer) timer.start();

        int nCols = seismicBoundingBox.nCols;
		float[] traceVals = new float[nSlices];
        for (int col = 0; col  < nCols; col++)
        {
			if (isCanceled()) return false;
            // getTraceData and output
            // traceVals = .....;
		    outputFloatBuffer.put(traceVals);
        }
        if (runTimer) timer.stopPrint("process " + nCols + " traces: ");

        return true;
    }

    public boolean isOutputDataFloat()
    {
		return seismicBoundingBox != null && inputVolumesAreFloat();
    }

	boolean inputVolumesAreFloat()
	{
//		for (int n = 0; n < inputVolumes.length; n++)
//			if (inputVolumes[n].rowFloatFilename == null)return false;
		return true;
	}

    public void finalizeBlockInput(int nBlock)
    {
        cleanInputMappedBuffers();
        super.finalizeBlockInput(nBlock);
    }

    public void finalizeVolumeInput()
    {
        clearInputMappedBuffers();
        cleanInputMappedBuffers();
        closeInputMappedBuffers();
    }

    public void finalizeVolumeOutput()
    {
    	super.finalizeVolumeOutput();
    	// remove histogram bias from spiked curvature values
    	outputVolume.dataHist[0] = 0.0f;
    }


    private void clearInputMappedBuffers()
    {
        if(outputFloatBuffer != null) outputFloatBuffer.clear0();
    }

    private void cleanInputMappedBuffers()
    {
        if(outputFloatBuffer != null) outputFloatBuffer.clean();
    }

    private void closeInputMappedBuffers()
    {
        if(outputFloatBuffer != null) outputFloatBuffer.close();
    }
}