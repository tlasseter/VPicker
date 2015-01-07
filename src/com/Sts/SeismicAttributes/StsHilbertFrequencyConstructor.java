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

public class StsHilbertFrequencyConstructor extends StsSeismicVolumeConstructor
{
    StsSeismicVolume data, hilbert;
    int nSamples;
    double[] traceDoubles, hilbertTraceDoubles;
    double[] phase;

    int windowHalfSize = 30;
    int windowSize = 2*windowHalfSize + 1;
    double maxFreq;

    StsTimer timer = null;
    boolean runTimer = false;

    static final double PI_OVER_2 = Math.PI/2;

    private StsHilbertFrequencyConstructor(StsModel model, StsSeismicVolume volume, boolean isDataFloat, StsProgressPanel panel)
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
            if (hilbert == null)return;
        }
        inputVolumes = new StsSeismicVolume[] {volume, hilbertVolume};
        volumeName = HILBERT_FREQ;
        String fullStemname = volume.stemname + "." + volumeName;
        maxFreq = 1000/volume.getZInc();
        outputVolume = StsSeismicVolume.initializeAttributeVolume(model, volume, 0.0f, (float)maxFreq, true, true, volume.stemname, volumeName, "rw");
//        initializeVolume(model, volume, 0.0f, (float)maxFreq, StsParameters.HILBERT_FREQ, volume.getZDomain());

		if(panel != null) panel.initialize(nInputRows * nInputCols);

        nSamples = volume.nSlices;
        traceDoubles = new double[nSamples];
        hilbertTraceDoubles = new double[nSamples];
        phase = new double[nSamples];
        createOutputVolume();
    }

    static public StsHilbertFrequencyConstructor constructor(StsModel model, StsSeismicVolume data, boolean isDataFloat, StsProgressPanel panel)
    {
        try
        {
            StsHilbertFrequencyConstructor hilbertFreqConstructor = new StsHilbertFrequencyConstructor(model, data, isDataFloat, panel);
            return hilbertFreqConstructor;
        }
        catch(Exception e)
        {
            StsMessage.printMessage("StsHilbertFrequencyConstructor.constructor() failed.");
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

    final private boolean processTrace(int nTrace, StsMappedFloatBuffer outputFloatBuffer, StsSeismicVolume vol)
    {
        int n = 0;

        try
        {
            double cycle, prevCycle = 0.0;
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
					double phase = Math.atan2(hilbert, data); // returns phase between -PI and +PI
                    phase += Math.PI; // move to range of 0 to 2PI
                    cycle = phase/(2*Math.PI); // convert to cycle range (0 to 1)
                    double dCycle;
                    int dCycleInt = 0;
                    if(n == 0)
                        dCycle = 0;
                    else
                    {
                        dCycle = cycle - prevCycle;
                        if(dCycle < 0) dCycle += 1;
                        dCycleInt = StsMath.minMax((int)(dCycle*254), 0, 254);
                    }
                    double freq = dCycle*maxFreq;
					outputFloatBuffer.put((float)freq);
				    vol.accumulateHistogram(dCycleInt);
				}
			}
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsHilbertFrequencyConstructor.processTrace() failed at  trace " + nTrace + " sample " + n,
                                         e, StsException.WARNING);
            return false;
        }
    }
}
