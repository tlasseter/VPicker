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
public  class StsCurvatureVolumeConstructor extends StsVolumeConstructor
{
    public StsPatchVolume patchVolume;
    boolean inputIsFloat = true;
	StsTimer timer = null;
	boolean runTimer = false;
	private int nSlices;
	String stsDirectory;
	String rowFloatFilename;
    int maxInterpolationSize;

    static final float nullValue = StsParameters.nullValue;

	//static methods
    static public StsSeismicVolume createInterpolatedVolume(StsModel model, StsPatchVolume inputVolume, StsProgressPanel panel, byte curvatureType, int maxInterpolationSize)
	{
    	int filterSize = inputVolume.filterSize;
    	String attributeName = StsSurfaceCurvatureAttribute.CURV_ATTRIBUTE_NAMES[curvatureType];
		if(panel != null) panel.appendLine("Creating " + attributeName + " volume.");

		StsSeismicVolume seismicVolume = checkGetAttributeVolume(model, inputVolume.getName() + "." + attributeName);
		if (seismicVolume != null)return seismicVolume;

		StsCurvatureVolumeConstructor attributeVolume = StsCurvatureVolumeConstructor.constructor(model, inputVolume, true, panel, curvatureType, filterSize, maxInterpolationSize);

		if(panel != null) panel.appendLine("Successfully created " + attributeName + " volume.");

		return attributeVolume.outputVolume;
	}

	static public StsSeismicVolume checkGetAttributeVolume(StsModel model, String attributeStemname)
	{
		StsSeismicVolume attributeVolume;
		try
		{
			attributeVolume = (StsSeismicVolume) model.getObjectWithName(StsSeismicVolume.class, attributeStemname);
			if (attributeVolume == null) return null;
			boolean deleteVolume = StsYesNoDialog.questionValue(model.win3d, "Volume " + attributeStemname + " already loaded. Delete and recreate?");
			if (!deleteVolume) return attributeVolume;
			attributeVolume.delete();
			attributeVolume = null;
			return null;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSeismicVolumeConstructor.checkVolumeExistence() failed.", e, StsException.WARNING);
			return null;
		}
	}

	static public StsSeismicVolume getExistingVolume(StsModel model, String attributeStemname)
	{
		try
		{
			return (StsSeismicVolume) model.getObjectWithName(StsSeismicVolume.class, attributeStemname);
		}
		catch (Exception e)
		{
			StsException.outputException("StsSeismicVolumeConstructor.checkVolumeExistence() failed.", e, StsException.WARNING);
			return null;
		}
	}


    static public StsCurvatureVolumeConstructor constructor(StsModel model, StsPatchVolume data,
    		boolean isDataFloat, StsProgressPanel panel, byte curvatureType, int filterSize, int maxInterpolationSize)
    {
        try
        {
            return new StsCurvatureVolumeConstructor(model, data, isDataFloat, panel, curvatureType, filterSize, maxInterpolationSize);
        }
        catch (Exception e)
        {
            StsMessage.printMessage("StsCurvatureVolumeConstructor.constructor() failed.");
            return null;
        }
    }

    //private constructor
    private StsCurvatureVolumeConstructor(StsModel model, StsPatchVolume volume,
    		boolean isDataFloat, StsProgressPanel panel, byte curvatureType, int filter, int maxInterpolationSize)
    {
        this.model = model;
        this.panel = panel;
        this.maxInterpolationSize = maxInterpolationSize;
        patchVolume = volume;

        StsSeismicBoundingBox seismicBoundingBox = new StsSeismicBoundingBox(patchVolume, false);
        //TODO simplify:  get from patchVolume.croppedBoundingBox
        // seismicBoundingBox.initializeToBoundingBox(volume.croppedBoundingBox);
    /*
        nInputRows = volume.nRows;
		nInputCols = volume.nCols;
        seismicBoundingBox.setRowNumMin(volume.getRowNumFromRow(volume.rowMin));
        seismicBoundingBox.setRowNumMax(volume.getRowNumFromRow(volume.rowMax));
        seismicBoundingBox.setColNumMin(volume.getColNumFromCol(volume.colMin));
        seismicBoundingBox.setColNumMax(volume.getColNumFromCol(volume.colMax));
        seismicBoundingBox.zMin = volume.zMin;
        seismicBoundingBox.zMax = volume.zMax;
    */
        seismicBoundingBox.zInc = patchVolume.interpolatedZInc;
        seismicBoundingBox.nSlices = patchVolume.nInterpolatedSlices;
        seismicBoundingBox.zDomain = patchVolume.zDomain;
        seismicBoundingBox.stsDirectory = patchVolume.stsDirectory;
        nSlices = patchVolume.nInterpolatedSlices;
        volumeName = StsSurfaceCurvatureAttribute.CURV_ATTRIBUTE_NAMES[(int)curvatureType]+filter+"x"+filter;
        float dataMin = patchVolume.getDataMin();
        float dataMax = patchVolume.getDataMax();
        outputVolume = StsSeismicVolume.initializeAttributeVolume(model, seismicBoundingBox, dataMin, dataMax, true, true, volume.getSeismicName(), volumeName, "rw");
        createOutputVolume();
    }


    public StsCurvatureVolumeConstructor()
	{
	}

	public StsCurvatureVolumeConstructor(StsPatchVolume inputVolume, String attributeName)
	{
        this.patchVolume = inputVolume;
		this.volumeName = attributeName;
	}

   public void initialize(String attributeName)
   {
       this.patchVolume = null;
       this.volumeName = attributeName;
   }


    public void createVolume(StsPatchVolume inputVolume)
	{
		this.patchVolume = inputVolume;
        createOutputVolume();
    }

    public void initializeVolumeInput()
    {
        if(patchVolume == null) return;
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
        return doProcessInputBlock(nBlock, null);
    }

    public boolean doProcessInputBlock(int nBlock, StsMappedBuffer[] inputBuffers)
    {
        if (runTimer) timer.start();

        int row = nBlock;
        int nCols = patchVolume.nCols;
		float[] traceVals = new float[nSlices];
        int[] patchRange = patchVolume.getPatchRangeForRow(row);
        for (int col = 0; col  < nCols; col++)
        {
			if (isCanceled()) return false;
			if (patchVolume.getTraceCurvature(row, col, traceVals, patchRange) )
			{
				interpolateTrace(traceVals, outputFloatBuffer, outputVolume, maxInterpolationSize);
			}
			else
			{
				outputFloatBuffer.put(traceVals);
			}
        }
        if (runTimer) timer.stopPrint("process " + nCols + " traces for block " + nBlock + ":");

        return true;
    }

    final private boolean interpolateTrace(float[] trace, StsMappedFloatBuffer floatBuffer, StsSeismicVolume vol, int maxInterpolationSize)
    {
    	if (maxInterpolationSize == 0)
    		maxInterpolationSize = nSlices /10;
        if(floatBuffer == null)return false;
        try
        {
            int prev = -1;
            float prevValue = nullValue;
            int next = 0;
            float nextValue;
            while(next != -1)
            {
                next = getNextFilledIndex(trace, nSlices, prev+1);
                if(next == -1)
                    nextValue = nullValue;
                else
                    nextValue = trace[next];
                fillIntermediateValues(trace, nSlices, prev, prevValue, next, nextValue, maxInterpolationSize);
                prev = next;
                prevValue = nextValue;
            }
            floatBuffer.put(trace);
            for (int n = 0; n < nSlices; n++)
				vol.accumulateHistogram((byte) trace[n]);
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsCurvatureVolumeConstructor.interpolateTrace() failed.",
                                         e, StsException.WARNING);
            return false;
        }
    }

    private int getNextFilledIndex(float[] trace, int nValues, int startIndex)
    {
           for (int n = startIndex; n < nValues; n++)
            	if (trace[n] != nullValue) return n;
           return -1;
    }

    private void fillIntermediateValues(float[] trace, int nValues, int prev, float prevValue, int next, float nextValue, int maxInterpolationSize)
    {
        if(prevValue == nullValue && nextValue == nullValue) return;
        boolean hasGap = next - prev > 2*maxInterpolationSize;
        if(prevValue == nullValue || hasGap)
        {
            prev = Math.max(0, next - maxInterpolationSize);
            for(int n = prev; n < next; n++)
                trace[n] = nextValue;
            return;
        }
        if(nextValue == nullValue || hasGap)
        {
            next = Math.min(nValues-1, prev + maxInterpolationSize);
            for(int n = prev+1; n < next; n++)
                trace[n] = prevValue;
            return;
        }
        // both values are not null and hasGap is false
        int nInc = next - prev;
        if(nInc == 0) return;
        float vInc = (nextValue - prevValue)/nInc;
        float v = prevValue + vInc;
        for(int n = prev+1; n < next; n++, v += vInc)
            trace[n] = v;
    }


    public boolean initializeBlockInput(int nBlock)
    {
        if (debug) System.out.println("Remapping inline file channel. outputPosition: " + outputPosition + " nOutputSamplesPerInputBlock: " + nOutputSamplesPerInputBlock);
        if(!outputFloatBuffer.map(outputPosition, nOutputSamplesPerInputBlock)) return false;
        outputPosition += nOutputSamplesPerInputBlock;

        return true;
    }

    public boolean isOutputDataFloat()
    {
		return patchVolume != null && inputVolumesAreFloat();
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
