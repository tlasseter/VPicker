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

public class StsHilbertAmplitudeConstructor extends StsSeismicVolumeConstructor
{
    StsSeismicVolume data, hilbert;
    int nSamples;
    double[] traceDoubles, hilbertTraceDoubles;

    int windowHalfSize = 30;
    int windowSize = 2 * windowHalfSize + 1;

    StsTimer timer = null;
    boolean runTimer = false;

    private StsHilbertAmplitudeConstructor(StsModel model, StsSeismicVolume volume, boolean isDataFloat, StsProgressPanel panel)
    {
        this.model = model;
        this.panel = panel;

		nInputRows = volume.nRows;
		nInputCols = volume.nCols;

        String hilbertStemname = new String(volume.stemname + "." + HILBERT);
        hilbert = (StsSeismicVolume) model.getObjectWithName(StsSeismicVolume.class,
            hilbertStemname);
        if (hilbert == null)
        {
            hilbert = StsSeismicVolume.checkLoadFromStemname(model, volume.stsDirectory, hilbertStemname, true);
        }
        if (hilbert == null)
        {
            StsHilbertTransformConstructor hilbertAttribute = StsHilbertTransformConstructor.constructor(model, volume, isDataFloat, panel);
            if (hilbertAttribute == null) return;
            hilbert = hilbertAttribute.getVolume();
        }

        inputVolumes = new StsSeismicVolume[] {volume, hilbert};
        volumeName = HILBERT_AMPLITUDE;
        outputVolume = StsSeismicVolume.initializeAttributeVolume(model, volume, volume.dataMin, volume.dataMax, true, true, volume.stemname, volumeName, "rw");

	    if(panel != null) panel.initialize(nInputRows * nInputCols);

        nSamples = volume.nSlices;
        traceDoubles = new double[nSamples];
        hilbertTraceDoubles = new double[nSamples];
        createOutputVolume();
    }

    static public StsHilbertAmplitudeConstructor constructor(StsModel model, StsSeismicVolume data, boolean isDataFloat, StsProgressPanel panel)
    {
        try
        {
            StsHilbertAmplitudeConstructor hilbertAmplitudeConstructor = new StsHilbertAmplitudeConstructor(model, data, isDataFloat, panel);
            return hilbertAmplitudeConstructor;
        }
        catch (Exception e)
        {
            StsMessage.printMessage("StsHilbertAmplitudeConstructor.constructor() failed.");
            return null;
        }
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
            inputBuffers[0].get(traceDoubles);
            inputBuffers[1].get(hilbertTraceDoubles);
            processTrace(nTrace, outputFloatBuffer, outputVolume);
        }
        if (runTimer)
        {
            timer.stopPrint("process " + nInputBlockTraces + " traces for block " + nBlock + ":");

        }
        return true;
    }

    final private boolean processTrace(int nTrace, StsMappedFloatBuffer outputFloatBuffer, StsSeismicVolume vol)
    {
        int n = 0;

        try
        {
            for (n = 0; n < nSamples; n++)
			{
				if (hilbertTraceDoubles[n] == StsParameters.nullValue || traceDoubles[n] == StsParameters.nullValue)
				{
					outputFloatBuffer.put(StsParameters.nullValue);
				}
				else
				{
					double data = traceDoubles[n];
					double hilbert = hilbertTraceDoubles[n];
					double result = Math.sqrt(data * data + hilbert * hilbert);
					outputFloatBuffer.put( (float) result);
					vol.accumulateHistogram((int)(2*result));
				}
			}
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsHilbertAmplitudeConstructor.processTrace() failed at  trace " + nTrace + " sample " +
                                         n,
                                         e, StsException.WARNING);
            return false;
        }
    }
}
