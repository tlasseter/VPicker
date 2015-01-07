package com.Sts.SeismicAttributes;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsVirtualVolumeMathConstructor extends StsVirtualVolumeConstructor
{
    StsSeismicVolume ouputVolume;
    int nSamples;
    float[] volOneDoubles, volTwoDoubles;

    StsTimer timer = null;
    boolean runTimer = false;

    private StsVirtualVolumeMathConstructor(StsModel model, StsSeismicVolume[] volumes, StsVirtualVolume vVolume, boolean isDataFloat, StsProgressPanel panel)
    {
        super(model, volumes, vVolume);
        this.panel = panel;

        nInputRows = volumes[0].nRows;
		nInputCols = volumes[0].nCols;
	    virtualVolume.setDataMax(-StsParameters.largeFloat);
	    virtualVolume.setDataMin(StsParameters.largeFloat);
	    
        //outputVolume = StsSeismicVolume.initializeAttributeVolume(model, vVolume, dataMin, dataMax, true, false, volumeName, "rw");
        outputVolume = StsSeismicVolume.initializeAttributeVolume(model, vVolume, vVolume.dataMin, vVolume.dataMax, true, false, volumeName, "rw");
	    if(panel != null) panel.initialize(nInputRows * nInputCols);	    
        createOutputVolume();
    }

    static public StsVirtualVolumeMathConstructor constructor(StsModel model, StsSeismicVolume[] data, StsVirtualVolume vVolume, boolean isDataFloat, StsProgressPanel panel)
    {
        try
        {
            StsVirtualVolumeMathConstructor mathVVConstructor = new StsVirtualVolumeMathConstructor(model, data, vVolume, isDataFloat, panel);
            return mathVVConstructor;
        }
        catch (Exception e)
        {
            StsMessage.printMessage("StsVirtualVolumeMathConstructor.constructor() failed.");
            return null;
        }
    }

    public boolean doProcessInputBlock(int nBlock, StsMappedBuffer[] inputBuffers)
    {
    	
    	if(volOneDoubles == null)
    	{
    		volOneDoubles = new float[nInputSamplesPerBlock];
    		if(virtualVolume.volumes.getSize() > 1)
    			volTwoDoubles = new float[nInputSamplesPerBlock];
    	}
        if (runTimer)
        {
            timer.start();
        }
		if (isCanceled()) return false;
		if(panel != null) 
			panel.setValue(nBlock);
		
        StsMappedBuffer inputBuffer = inputBuffers[0];
        inputBuffer.position(0);
        inputBuffers[0].get(volOneDoubles);
        float[] result;
		if(inputBuffers.length > 1)		
		{
			inputBuffers[1].get(volTwoDoubles);	
			result =  ((StsMathVirtualVolume)virtualVolume).processBlock(volOneDoubles, volTwoDoubles, nInputSamplesPerBlock);
		}
		else
			result =  ((StsMathVirtualVolume)virtualVolume).processBlock(volOneDoubles, nInputSamplesPerBlock);						
	    outputFloatBuffer.put(result);

	    outputVolume.setDataMax(virtualVolume.dataMax);
	    outputVolume.setDataMin(virtualVolume.dataMin);
	    outputVolume.writeHeaderFile();
	    outputVolume.reinitializeColorscale();
        if (runTimer)
        {
            timer.stopPrint("process " + nInputBlockTraces + " traces for block " + nBlock + ":");

        }
        return true;
    }
}
