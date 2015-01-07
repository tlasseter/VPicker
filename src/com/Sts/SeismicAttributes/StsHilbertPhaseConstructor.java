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

public class StsHilbertPhaseConstructor extends StsSeismicVolumeConstructor
{
    StsSeismicVolume data, hilbert;
    int nSamples;
    double[] traceDoubles, hilbertTraceDoubles;

    int windowHalfSize = 30;
    int windowSize = 2 * windowHalfSize + 1;

    StsTimer timer = null;
    boolean runTimer = false;

    private StsHilbertPhaseConstructor(StsModel model, StsSeismicVolume volume, boolean isDataFloat, StsProgressPanel panel)
    {
        this.model = model;
        this.panel = panel;

		nInputRows = volume.nRows;
		nInputCols = volume.nCols;

        String hilbertStemname = new String(volume.stemname + "." + HILBERT);
        StsSeismicVolume hilbertVolume = (StsSeismicVolume) model.getObjectWithName(StsSeismicVolume.class, hilbertStemname);
        if (hilbertVolume == null)
        {
            hilbertVolume = StsSeismicVolume.checkLoadFromStemname(model, volume.stsDirectory, hilbertStemname, true);
        }
        if (hilbertVolume == null)
        {
            StsHilbertTransformConstructor hilbert = StsHilbertTransformConstructor.constructor(model, volume, isDataFloat, panel);
            if (hilbert == null) return;
        }

        inputVolumes = new StsSeismicVolume[] {volume, hilbertVolume};
        volumeName = HILBERT_PHASE;
        outputVolume = StsSeismicVolume.initializeAttributeVolume(model, volume, -(float) Math.PI, (float) Math.PI, true, true, volume.stemname, volumeName, "rw");
//        initializeVolume(model, volume, -(float) Math.PI, (float) Math.PI, StsParameters.HILBERT_PHASE, volume.getZDomain());

		if(panel != null) panel.initialize(nInputRows * nInputCols);

        nSamples = volume.nSlices;
        traceDoubles = new double[nSamples];
        hilbertTraceDoubles = new double[nSamples];
        createOutputVolume();
    }

    static public StsHilbertPhaseConstructor constructor(StsModel model, StsSeismicVolume data, boolean isDataFloat, StsProgressPanel panel)
    {
        try
        {
            return new StsHilbertPhaseConstructor(model, data, isDataFloat, panel);
        }
        catch (Exception e)
        {
            StsMessage.printMessage("StsHilbertPhaseConstructor.constructor() failed.");
            return null;
        }
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

    final private boolean processTrace(int nTrace, StsMappedFloatBuffer floatBuffer, StsSeismicVolume vol)
    {
        int n = 0;

        try
        {
            for (n = 0; n < nSamples; n++)
            {
				if(hilbertTraceDoubles[n] == StsParameters.nullValue || traceDoubles[n] == StsParameters.nullValue)
				{
					floatBuffer.put(StsParameters.nullValue);
				}
				else
				{
					double phase = Math.atan2(hilbertTraceDoubles[n], traceDoubles[n]);
					if(floatBuffer != null) floatBuffer.put( (float) phase);
					int phaseInt = (int) (phase * 127 / Math.PI + 127);
					vol.accumulateHistogram(phaseInt);
				}
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsHilbertPhaseConstructor.processTrace() failed at  trace " + nTrace + " sample " + n,
                                         e, StsException.WARNING);
            return false;
        }
    }
}
