package com.Sts.Utilities.Seismic;

import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.Utilities.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Aug 5, 2009
 * Time: 9:55:34 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsHilbertTransform
{
    static double[] hilbertFilter = null;
    static final int defaultWindowHalfSize = 30;
    static int windowHalfSize = defaultWindowHalfSize;
    static int windowSize = 2 * windowHalfSize + 1;

    StsTimer timer = null;
    boolean runTimer = false;

    static public void staticInitialize()
    {
        staticInitialize(defaultWindowHalfSize);
    }

    static public void staticInitialize(int newWindowHalfSize)
    {
        if(windowHalfSize != newWindowHalfSize)
        {
            windowHalfSize = newWindowHalfSize;
            windowSize = 2 * windowHalfSize + 1;
        }

        if(hilbertFilter != null && hilbertFilter.length == newWindowHalfSize) return;

        hilbertFilter = new double[windowSize];
        hilbertFilter[windowHalfSize] = 0;
        for (int i = 1; i <= windowHalfSize; i++)
        {
            double taper = 0.54 + 0.46 * Math.cos(Math.PI * (double) i / (double) (windowHalfSize));
            hilbertFilter[windowHalfSize + i] = (taper * ( - (double) (i % 2) * 2.0 / (Math.PI * (double) (i))));
            hilbertFilter[windowHalfSize - i] = -hilbertFilter[windowHalfSize + i];
        }
    }

    static public double[] computeHilbert(double[] traceDoubles)
    {
        int nSamples = traceDoubles.length;
        double[] hilbert = new double[nSamples];
        StsConvolve.convolve(traceDoubles, hilbertFilter, hilbert, nSamples, 0, nSamples-1, windowSize, windowHalfSize);
        return hilbert;
    }

    static public final float[] computeHilbert(float[] traceFloats)
    {
        int nSamples = traceFloats.length;
        float[] hilbert = new float[nSamples];
        StsConvolve.convolve(traceFloats, hilbertFilter, hilbert, nSamples, 0, nSamples-1, windowSize, windowHalfSize);
        return hilbert;
    }

    static public final double[] computePhase(double[] hilbert, double[] data)
    {
        int nSamples = data.length;
        double[] phase = new double[nSamples];
        for(int n = 0; n < nSamples; n++)
            phase[n] = Math.atan2(hilbert[n], data[n]);
        return phase;
    }

    static public final float[] computePhase(float[] hilbert, float[] data)
    {
        int nSamples = data.length;
        float[] phase = new float[nSamples];
        for(int n = 0; n < nSamples; n++)
            phase[n] = (float)Math.atan2(hilbert[n], data[n]);
        return phase;
    }

    static public final double[] computeAmplitude(double[] hilbert, double[] data)
    {
        int nSamples = data.length;
        double[] amplitude = new double[nSamples];
        for(int n = 0; n < nSamples; n++)
        {
            double d = data[n];
			double h = hilbert[n];
			amplitude[n] = Math.sqrt(d * d + h * h);
        }
        return amplitude;
    }

    static public final float[] computeAmplitude(float[] hilbert, float[] data)
    {
        int nSamples = data.length;
        float[] amplitude = new float[nSamples];
        for(int n = 0; n < nSamples; n++)
        {
            double d = data[n];
			double h = hilbert[n];
			amplitude[n] = (float)Math.sqrt(d * d + h * h);
        }
        return amplitude;
    }
}
