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

public class StsHilbertTransformConstructor extends StsSeismicVolumeConstructor
{
	StsSeismicVolume inputVolume;
	int nSamples;
	double[] hilbertFilter = null;
	double[] transformValues;
	double[] traceDoubles = null;

	double scaleFactor, scaleOffset;

	static int windowHalfSize = 30;
	static int windowSize = 2 * windowHalfSize + 1;

	StsTimer timer = null;
	boolean runTimer = false;

	private StsHilbertTransformConstructor(StsModel model, StsSeismicVolume inputVolume, boolean isDataFloat, StsProgressPanel panel)
	{
		super(new StsSeismicVolume[] {inputVolume}, HILBERT);
        this.model = model;
        this.panel = panel;
		if (runTimer) timer = new StsTimer();
        volumeName = HILBERT;
        outputVolume = StsSeismicVolume.initializeAttributeVolume(model, inputVolume, inputVolume.dataMin, inputVolume.dataMax, true, true, inputVolume.stemname, volumeName, "rw");
//        initializeVolume(model, inputVolume, inputVolume.dataMin, inputVolume.dataMax, StsParameters.HILBERT_TRANS, inputVolume.getZDomain());
		this.inputVolume = inputVolume;
		nSamples = inputVolume.nSlices;
		nInputRows = inputVolume.nRows;
		nInputCols = inputVolume.nCols;
		transformValues = new double[nSamples];
		traceDoubles = new double[nSamples];

		hilbertFilter = new double[windowSize];
		hilbertFilter[windowHalfSize] = 0;
		for (int i = 1; i <= windowHalfSize; i++)
		{
			double taper = 0.54 + 0.46 * Math.cos(Math.PI * (double) i / (double) (windowHalfSize));
			hilbertFilter[windowHalfSize + i] = (taper * ( - (double) (i % 2) * 2.0 / (Math.PI * (double) (i))));
			hilbertFilter[windowHalfSize - i] = -hilbertFilter[windowHalfSize + i];
		}
		outputVolume.initializeScaling();
		if(panel != null) panel.initialize(nInputRows * nInputCols);
        createOutputVolume();
    }

   static public StsHilbertTransformConstructor constructor(StsModel model, StsSeismicVolume inputVolume, boolean isDataFloat, StsProgressPanel panel)
    {
        try
        {
            return new StsHilbertTransformConstructor(model, inputVolume, isDataFloat, panel);
        }
        catch(Exception e)
        {
            StsMessage.printMessage("StsHilbertTransformConstructor.constructor() failed.");
            return null;
        }
    }

	private void initializeScaling()
	{
		scaleFactor = 254/(outputVolume.dataMax - outputVolume.dataMin);
		scaleOffset = -outputVolume.dataMin*scaleFactor;
	}

    public boolean doProcessInputBlock(int nBlock, StsMappedBuffer[] inputBuffers)
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
            inputBuffer.get(traceDoubles);
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
            StsConvolve.convolve(traceDoubles, hilbertFilter, transformValues, nSamples, 0, nSamples-1, windowSize, windowHalfSize);

            for (n = 0; n < nSamples; n++)
            {
				int scaledValue = (int)(transformValues[n]*scaleFactor + scaleOffset);
                scaledValue = StsMath.minMax(scaledValue, 0, 254);
                volume.accumulateHistogram(scaledValue);
                byte b = StsMath.unsignedIntToUnsignedByte(scaledValue);
				floatBuffer.put((float)transformValues[n]);
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
