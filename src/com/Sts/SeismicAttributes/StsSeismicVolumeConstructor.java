package com.Sts.SeismicAttributes;

import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.util.*;

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
public abstract class StsSeismicVolumeConstructor extends StsVolumeConstructor
{
	/** random access file for reading input segy file */
//    transient private RandomAccessFile inputRowFile;
//    transient private FileChannel inputRowChannel;
//    transient private MappedByteBuffer inputBuffer = null;

    public StsSeismicVolume[] inputVolumes;
    StsMappedBuffer[] inputBuffers = null;
    boolean inputIsFloat = false;

    static public final String HILBERT = "hilbert";
	static public final String HILBERT_PHASE = "hilbertPhase";
	static public final String HILBERT_AMPLITUDE = "hilbertAmplitude";
	static public final String HILBERT_FREQ = "hilbertFreq";
    static public final String FREQ_ENHANCE = "freqEnhance";
    static public final String ANALOGUE = "analog";
    static public final String INTERPOLATE = "interpolate";

    public  boolean doProcessInputBlock(int nBlock, StsMappedBuffer[] inputBuffers, String mode) { return doProcessInputBlock(nBlock, inputBuffers);}

	abstract public  boolean doProcessInputBlock(int nBlock, StsMappedBuffer[] inputBuffers);

    public StsSeismicVolumeConstructor()
	{
	}

	public StsSeismicVolumeConstructor(StsSeismicVolume[] inputVolumes, String attributeName)
	{
        this.inputVolumes = inputVolumes;
		this.volumeName = attributeName;
	}

    public void initialize(StsSeismicVolume inputVolume, String attributeName)
    {
        this.inputVolumes = new StsSeismicVolume[] { inputVolume };
        this.volumeName = attributeName;
   }

   public void initialize(String attributeName)
   {
       this.inputVolumes = null;
       this.volumeName = attributeName;
   }

    static public String[] getPossibleAttributeNames() { return new String[] { HILBERT, HILBERT_PHASE, HILBERT_AMPLITUDE, HILBERT_FREQ}; }

	static public String[] getAvailableAttributeNames(StsSeismicVolume volume, StsModel model)
	{
		String[] attributeNames = getPossibleAttributeNames();
		int nAttributes = attributeNames.length;
		String[] availableNames = new String[nAttributes];
		int nAvailable = 0;
		for (int n = 0; n < nAttributes; n++)
		{
			String name = new String(volume.stemname + "." + attributeNames[n]);
			StsSeismicVolume attributeVolume = (StsSeismicVolume) model.getObjectWithName(StsSeismicVolume.class, name);
			if (attributeVolume != null) availableNames[nAvailable++] = attributeNames[n];
		}
		return (String[]) StsMath.trimArray(availableNames, nAvailable);
	}

	static public void createAttributeVolumes(StsModel model, StsSeismicVolume inputVolume, ArrayList attributeNames, boolean isDataFloat, StsProgressPanel panel)
	{
		for (int n = 0; n < attributeNames.size(); n++)
		{
			if(panel != null) panel.appendLine("Creating " + (String)attributeNames.get(n) + " volume.");
//			panel.appendLine(new Date(System.currentTimeMillis()).toString() + "\n     Creating " +
//							 (String) attributeNames.get(n) + " volume.");
			createAttributeVolume(model, inputVolume, (String) attributeNames.get(n), isDataFloat, panel);
			if(panel != null) panel.appendLine("Successfully created " + (String)attributeNames.get(n) + " volume.");
//			panel.appendLine(new Date(System.currentTimeMillis()).toString() + "\n     Successfully created " +
//							 (String) attributeNames.get(n) + " volume.\n");
		}
	}

    static public StsSeismicVolume createAttributeVolume(StsModel model, StsSeismicVolume inputVolume,
		String attributeName, boolean isDataFloat, StsProgressPanel panel)
	{
		StsSeismicVolumeConstructor attributeVolume = null;
		StsSeismicVolume seismicVolume;

		seismicVolume = checkGetAttributeVolume(model, inputVolume.stemname + "." + attributeName);
		if (seismicVolume != null)return seismicVolume;

		if (attributeName.equals(HILBERT))
			attributeVolume = StsHilbertTransformConstructor.constructor(model, inputVolume, isDataFloat, panel);
		else if (attributeName.equals(HILBERT_PHASE))
			attributeVolume = StsHilbertPhaseConstructor.constructor(model, inputVolume, isDataFloat, panel);
		else if (attributeName.equals(HILBERT_AMPLITUDE))
			attributeVolume = StsHilbertAmplitudeConstructor.constructor(model, inputVolume, isDataFloat, panel);
		else if (attributeName.equals(HILBERT_FREQ))
			attributeVolume = StsHilbertFrequencyConstructor.constructor(model, inputVolume, isDataFloat, panel);
		if (attributeVolume == null)return null;

		seismicVolume = attributeVolume.getVolume();
		if(!seismicVolume.initialize(model))
		{
			seismicVolume.delete();
			return null;
		}
		return seismicVolume;
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

    public void createVolume(StsSeismicVolume[] inputVolumes)
	{
		this.inputVolumes = inputVolumes;
        createOutputVolume();
    }

    public void initializeVolumeInput()
    {
        super.initializeVolumeInput();
        outputFloatBuffer = outputVolume.createMappedFloatRowBuffer();
        outputPosition = 0;
        if(inputVolumes == null) return;
        int nInputVolumes = inputVolumes.length;
        this.inputIsFloat = inputVolumesAreFloat();
        if (nInputVolumes > 0)
        {
            inputBuffers = new StsMappedBuffer[nInputVolumes];
            for (int n = 0; n < nInputVolumes; n++)
                inputBuffers[n] = inputVolumes[n].createMappedRowBuffer("r");
        }
        inputPosition = 0;
    }

    public boolean allocateMemory()
    {
        int nInputVolumes = 0;
        if(inputVolumes != null)
            nInputVolumes = inputVolumes.length;
        if(nInputVolumes == 0)
            memoryAllocation = StsMemAllocVolumeProcess.constructor(outputVolume);
        else
            memoryAllocation = StsMemAllocVolumeProcess.constructor(inputVolumes[0], outputVolume);
        if(memoryAllocation == null) return false;
        nOutputRowsPerBlock = memoryAllocation.nOutputRowsPerBlock;
        nOutputBlocks = memoryAllocation.nOutputBlocks;
        nOutputSamplesPerInputBlock = memoryAllocation.nOutputSamplesPerInputBlock;
        nOutputSamplesPerRow = memoryAllocation.nOutputSamplesPerRow;
        return true;
    }

    public  boolean processBlockInput(int nBlock, String mode)
    {
        return doProcessInputBlock(nBlock, inputBuffers, mode);
    }

	public  boolean processBlockInput(int nBlock)
	  {
		  return doProcessInputBlock(nBlock, inputBuffers);
	  }

    public boolean initializeBlockInput(int nBlock)
    {
        if(!super.initializeBlockInput(nBlock)) return false;

        if(inputVolumes != null)
        {
            // System.out.println("Remapping inputBuffers. inputPosition: " + inputPosition + " nBlockInputSamples: " + nInputSamplesPerBlock);
            for (int n = 0; n < inputVolumes.length; n++)
                if(!inputBuffers[n].map(inputPosition, nInputSamplesPerBlock)) return false;
            inputPosition += nInputSamplesPerBlock;
        }

        if (debug) System.out.println("Remapping inline file channel. outputPosition: " + outputPosition + " nOutputSamplesPerInputBlock: " + nOutputSamplesPerInputBlock);
        if(!outputFloatBuffer.map(outputPosition, nOutputSamplesPerInputBlock)) return false;
        outputPosition += nOutputSamplesPerInputBlock;

        return true;
    }

    public boolean isOutputDataFloat()
    {
		return inputVolumes != null && inputVolumesAreFloat();
    }

	boolean inputVolumesAreFloat()
	{
		for (int n = 0; n < inputVolumes.length; n++)
			if (inputVolumes[n].rowFloatFilename == null)return false;
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

    private void clearInputMappedBuffers()
    {
        if (inputBuffers != null)
		{
			for (int n = 0; n < inputVolumes.length; n++)
				if (inputBuffers[n] != null) inputBuffers[n].clear0();
		}
        if(outputFloatBuffer != null) outputFloatBuffer.clear0();
    }

    private void cleanInputMappedBuffers()
    {
        if (inputBuffers != null)
		{
			for (int n = 0; n < inputVolumes.length; n++)
				if (inputBuffers[n] != null) inputBuffers[n].clean();
		}
        if(outputFloatBuffer != null) outputFloatBuffer.clean();
    }

    private void closeInputMappedBuffers()
    {
        if (inputBuffers != null)
		{
			for (int n = 0; n < inputVolumes.length; n++)
				if (inputBuffers[n] != null) inputBuffers[n].close();
		}
        if(outputFloatBuffer != null) outputFloatBuffer.close();
    }
}
