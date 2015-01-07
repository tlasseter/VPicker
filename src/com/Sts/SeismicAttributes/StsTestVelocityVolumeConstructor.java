package com.Sts.SeismicAttributes;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

import com.Sts.Actions.Wizards.Velocity.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsTestVelocityVolumeConstructor extends StsSeismicVolumeConstructor
{
    int nSamples;
    double[] traceDoubles;
	double[] velocities;
    float dataMin = 4.0f;
    float dataMax = 5.0f;
    double scaleFactor, scaleOffset;
    double tMin, tMax, tInc;

    StsTimer timer = null;
    boolean runTimer = false;
    public StsEditVelocityPanel editVelocityPanel;

//    StsBooleanLock booleanLock;

    static public final String TEST_VELOCITY = "testVelocity";

    private StsTestVelocityVolumeConstructor(StsModel model, StsSeismicVolume volume, boolean isDataFloat,
                                StsEditVelocityPanel editVelocityPanel)
    {
        super(new StsSeismicVolume[] {volume}, TEST_VELOCITY);
        this.model = model;
        tMin = volume.zMin;
        tMax = volume.zMax;
        tInc = volume.zInc;
        this.editVelocityPanel = editVelocityPanel;
//        this.panel = editVelocityPanel.progressPanel;
        if (runTimer) timer = new StsTimer();
        outputVolume = new StsSeismicVolume(false);
        outputVolume.initializeVolume(model, volume, dataMin, dataMax, true, true, TEST_VELOCITY, StsProject.TD_TIME, true, "rw");
        nSamples = volume.nSlices;
        nInputRows = volume.nRows;
        nInputCols = volume.nCols;
        panel.initialize(nInputRows * nInputCols);
        traceDoubles = new double[nSamples];
		velocities = new double[nSamples];
        initializeScaling();
        createOutputVolume();
    }

    public void createOutputVolume()
    {
        Runnable runnable = new Runnable()
        {
            public void run()
            {
                runCreateOutputVolume();
            }
        };
        StsToolkit.runRunnable(runnable);
    }

    protected boolean runCreateOutputVolume()
    {
        super.runCreateOutputVolume();
        StsTestVelocityVolumeConstructor.this.editVelocityPanel.setVolume(outputVolume);
        panel.finished();
        return true;
    }

    static public StsTestVelocityVolumeConstructor constructVolume(StsModel model, StsSeismicVolume inputVolume, boolean isDataFloat, StsEditVelocityPanel editVelocityPanel)
     {
         try
         {
             return new StsTestVelocityVolumeConstructor(model, inputVolume, isDataFloat, editVelocityPanel);
         }
         catch(Exception e)
         {
             StsMessage.printMessage("StsTestVelocityVolumeConstructor.constructor() failed.");
             return null;
         }
    }
    private void initializeScaling()
    {
        scaleFactor = 1.0/(inputVolumes[0].dataMax - inputVolumes[0].dataMin);
        scaleOffset = -inputVolumes[0].dataMin*scaleFactor;
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
            if(panel != null && (++nTracesDone%1000) == 0)
            {
                panel.setValue(nTracesDone);
                panel.appendLine("Processed " + nTracesDone);
            }
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
		double sum = 0;
		try
		{
			double t = tMin;
			double value = (traceDoubles[0]*scaleFactor + scaleOffset); // should be between 0 and 1
			double dataRange = dataMax - dataMin;
			double variableFraction = 0.2*dataRange;

			double rmsVelocity = dataMin;
			putVelocity(rmsVelocity, volume, floatBuffer);
			sum = rmsVelocity*rmsVelocity*t;
			double nextdSum = rmsVelocity*rmsVelocity*tInc/2;
			for (n = 1; n < nSamples; n++)
			{
				t += tInc;
				double prevdSum = nextdSum;
				value = (float)(traceDoubles[n]*scaleFactor + scaleOffset); // should be between 0 and 1
				double f = (t - tMin)/(tMax - tMin);

				double velocity = dataMin + f*dataRange*(0.8 + 0.2*value);
				nextdSum = velocity*velocity*tInc/2;
				sum += (prevdSum + nextdSum);
				rmsVelocity = Math.sqrt(sum/t);
				putVelocity(rmsVelocity, volume, floatBuffer);
			}
			return true;
		}
		catch (Exception e)
		{
			StsException.outputException("StsTestVelocityVolumeConstructor.processTrace() failed at  trace " + nTrace + " sample " + n,
										 e, StsException.WARNING);
			return false;
		}
    }
/*
	test version
    final private boolean processTrace(int nTrace, StsMappedByteBuffer inlineBuffer, StsMappedFloatBuffer inlineFloatBuffer, StsSeismicVolume volume)
    {
		int n = 0;
        try
        {
			double t = tMin;
			double value = (traceDoubles[0] * scaleFactor + scaleOffset); // should be between 0 and 1
			double dataRange = dataMax - dataMin;
			double variableFraction = 0.2 * dataRange;
			double velocity = dataMin;
			double sum = velocity;
			velocities[n] = velocity;
//			putVelocity(velocity, volume, inlineFloatBuffer);
			for (n = 1; n < nSamples; n++)
			{
				t += tInc;
				value = (float) (traceDoubles[n] * scaleFactor + scaleOffset); // should be between 0 and 1
				double f = (t - tMin) / (tMax - tMin);
				velocity = dataMin + dataRange * (f + 0.3 * (value - 0.5));
				sum += velocity;
				float avgVelocity = (float) (sum / (n + 1));
//				if(n == 200) System.out.println("vel " + velocity + " avg vel " + avgVelocity);
//				putVelocity(velocity, volume, inlineFloatBuffer);
				velocities[n] = avgVelocity;
			}
			double v1 = velocities[0];
			putVelocity(v1, volume, inlineFloatBuffer);
			for(n = 1; n < nSamples; n++)
			{
				double v0 = v1;
				v1 = velocities[n];
				double vi = v1*n - v0*(n-1);
				putVelocity(vi, volume, inlineFloatBuffer);
			}
/*
            double rmsVelocity = dataMin;
            putVelocity(rmsVelocity, volume, inlineFloatBuffer);
            sum = rmsVelocity*rmsVelocity*t;
            double nextdSum = rmsVelocity*rmsVelocity*tInc/2;
            for (n = 1; n < nSamples; n++)
            {
                t += tInc;
                double prevdSum = nextdSum;
                value = (float)(traceDoubles[n]*scaleFactor + scaleOffset); // should be between 0 and 1
                double f = (t - tMin)/(tMax - tMin);

                double velocity = dataMin + f*dataRange*(0.5 + 0.5*value);
                nextdSum = velocity*velocity*tInc/2;
                sum += (prevdSum + nextdSum);
                rmsVelocity = Math.sqrt(sum/t);
                putVelocity(rmsVelocity, volume, inlineFloatBuffer);
            }
*/
/*
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsTestVelocityVolumeConstructor.processTrace() failed at  trace " + nTrace + " sample " + n,
                                         e, StsException.WARNING);
            return false;
        }
    }
*/
    private void putVelocity(double velocity, StsSeismicVolume volume, StsMappedFloatBuffer floatBuffer)
    {
        double value = velocity - dataMin;
        int scaledValue = (int)(value*254);
        scaledValue = StsMath.minMax(scaledValue, 0, 254);
        volume.accumulateHistogram(scaledValue);
        byte b = StsMath.unsignedIntToUnsignedByte(scaledValue);
        if(floatBuffer != null) floatBuffer.put((float)velocity);
    }
}
