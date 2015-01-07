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

abstract public class StsVolumeInterpolationConstructor extends StsSeismicVolumeConstructor
{
	StsSeismicVolume inputVolume;
	int nSamples;

	double scaleFactor, scaleOffset;

	StsTimer timer = null;
	boolean runTimer = false;

	private StsVolumeInterpolationConstructor(StsModel model, StsSeismicVolume inputVolume, boolean isDataFloat, StsProgressPanel panel)
	{
		super(new StsSeismicVolume[] {inputVolume}, INTERPOLATE);
        this.model = model;
        this.panel = panel;
		if (runTimer) timer = new StsTimer();
        volumeName = INTERPOLATE;
        outputVolume = StsSeismicVolume.initializeAttributeVolume(model, inputVolume, inputVolume.dataMin, inputVolume.dataMax, true, true, inputVolume.stemname, volumeName, "rw");
		this.inputVolume = inputVolume;
		nSamples = inputVolume.nSlices;
		nInputRows = inputVolume.nRows;
		nInputCols = inputVolume.nCols;
		outputVolume.initializeScaling();
		if(panel != null) panel.initialize(nInputRows * nInputCols);
        createOutputVolume();
    }

   static public StsVolumeInterpolationConstructor constructor(StsModel model, StsSeismicVolume inputVolume, boolean isDataFloat, StsProgressPanel panel)
    {
        /*
        try
        {

            //return new StsVolumeInterpolationConstructor(model, inputVolume, isDataFloat, panel);
            return null;


        }
        catch(Exception e)
        */

        {
            StsMessage.printMessage("StsVolumeInterpolationConstructor.constructor() failed.");
            return null;
        }

    }

	private void initializeScaling()
	{
		scaleFactor = 254/(outputVolume.dataMax - outputVolume.dataMin);
		scaleOffset = -outputVolume.dataMin*scaleFactor;
	}

    public boolean doProcessInputBlock(int nBlock, StsMappedBuffer[] inputBuffers
    )
    {
        int nTrace = -1;

        if (runTimer)
        {
            timer.start();

        }
        StsMappedBuffer inputBuffer = inputBuffers[0];
        inputBuffer.position(0);
        for (nTrace = 0; nTrace < nInputBlockTraces; nTrace++)
        {
			if (isCanceled()) return false;
			if(panel != null && (++nTracesDone%1000) == 0) panel.setValue(nTracesDone);
//            inputBuffer.get(traceDoubles);
            processTrace(nTrace, outputFloatBuffer, outputVolume);
            if(debug && nTrace%1000 == 0) System.out.println("    processed trace: " + nTrace);
        }
        if (runTimer)
        {
            timer.stopPrint("process " + nInputBlockTraces + " traces for block " + nBlock + ":");

        }
        return true;
    }

    final private boolean processTrace(int nTrace, StsMappedFloatBuffer floatBuffer, StsSeismicVolume volume)
    {
        int n = 0;

        try
        {
 //           StsConvolve.convolve(traceDoubles, hilbertFilter, transformValues, nSamples, 0, nSamples-1, windowSize, windowHalfSize);

            for (n = 0; n < nSamples; n++)
            {
				int scaledValue = 0;
                // (int)(transformValues[n]*scaleFactor + scaleOffset);
                scaledValue = StsMath.minMax(scaledValue, 0, 254);
                volume.accumulateHistogram(scaledValue);
                byte b = StsMath.unsignedIntToUnsignedByte(scaledValue);
				// floatBuffer.put((float)transformValues[n]);
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsHilbertTransformConstructor.processTrace() failed at  trace " + nTrace + " sample " + n,
                                         e, StsException.WARNING);
            return false;
        }
    }
}