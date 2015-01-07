package com.Sts.Utilities.Filter;

import com.Sts.Utilities.*;

public class AGC
{

    /**
     * AGC - Automatic Gain Control
     * <p/>
     * Method to take full gather as input and then pass to version that takes trace as input (iterating over traces).
     *
     * @param gather      Input gather array (float array of trace samples)
     * @param nFirstTrace Index of trace to begin AGC application (starting with 0)
     * @param nLastTrace  Index of trace to end AGC application (starting with 0)
     * @param traceLength Number of samples in each trace contained within "gather"
     * @param windowWidth Number of samples for sliding AGC window
     * @param maxAmp      This is a bit confusing - it really should be "desired RMS amplitude for output" (samples will exceed this value)
     * @return Output is gather array with AGC applied
     */
	
	static boolean agcDebug = false;
	static boolean fastAGC = true;
	
    public static float[] instantaneousAGC(float[] gather, int nFirstTrace, int nLastTrace, int traceLength, int windowWidth, float maxAmp)
    {

        if (gather == null) return null;

        if ((nLastTrace - nFirstTrace) * traceLength > gather.length)
        {
            System.err.println("AGC: Number of traces greater than gather length - first trace: " + nFirstTrace +
                " last trace: " + nLastTrace + " gather traces: " + gather.length / traceLength);
            return gather;
        }

        float[] trace = new float[traceLength];         //current trace being AGC'ed
        float[] outGather = new float[gather.length];   //gather with AGC applied
        int position;                                   //current sample position in gather

        for (int nTrace = nFirstTrace; nTrace <= nLastTrace; nTrace++)
        {
            position = nTrace * traceLength;
            System.arraycopy(gather, position, trace, 0, traceLength);  //grab current trace from gather
            trace = instantaneousAGC(trace, windowWidth, maxAmp);   //apply agc
            if (trace != null)
            {
                System.arraycopy(trace, 0, outGather, position, traceLength); //put trace into output gather
            }

        }
        return outGather;
    }

    public static float[] rmsAGC(float[] gather, int nFirstTrace, int nLastTrace, int traceLength, int halfWindow, float maxAmp)
    {

        if (gather == null) return null;

        if ((nLastTrace - nFirstTrace) * traceLength > gather.length)
        {
            System.err.println("AGC: Number of traces greater than gather length - first trace: " + nFirstTrace +
                " last trace: " + nLastTrace + " gather traces: " + gather.length / traceLength);
            return gather;
        }

        float[] trace = new float[traceLength];         //current trace being AGC'ed
        float[] outGather = new float[gather.length];   //gather with AGC applied
        int position;                                   //current sample position in gather

        for (int nTrace = nFirstTrace; nTrace <= nLastTrace; nTrace++)
        {
            position = nTrace * traceLength;
            System.arraycopy(gather, position, trace, 0, traceLength);  //grab current trace from gather
            trace = rmsAGC(trace, halfWindow, maxAmp);   //apply agc
            if (trace != null)
            {
                System.arraycopy(trace, 0, outGather, position, traceLength); //put trace into output gather
            }

        }
        return outGather; 
    }

    /**
     * AGC - Automatic Gain Control
     * Algorithm based on Yilmaz, Seismic Data Processing, 1987 page 61
     * Implementation mimics behavior of Paradigm's Focus AGC (handling of boundary conditions, etc.)
     * <p/>
     * Summary:
     * Calculates average absolute amplitude for a given window.
     * The gain for that window will be = sFactor/avg_abs_amp.
     * The routine then moves the window up one sample and runs the calculation again.
     * <p/>
     * Once all of the gains are computed, data samples are multiplied by the gain value for each sample.
     * Samples at the beginning and end of the trace are given the gain from the nearest sample that
     * had a full window (i.e. - no partial windows are used).
     * <p/>
     * Error Handling:
     * If zero sFactor encountered, uses internally calculated sFactor (sets to reasonable value)
     * If window is greater than trace length, returns data unchanged.
     *
     * @param trace   Input trace samples (float array)
     * @param window  Length of sliding window (in samples)
     * @param sFactor Normalization amplitude (what you want the RMS amplitude of the output to be)
     * @return Output trace samples with AGC applied
     */
    public static float[] slowInstantaneousAGC(float[] trace, int window, double sFactor)
    {
        int length = trace.length;                    //length of input data array
        double[] scalars = new double[length];        //array off scalars to apply to data
        float[] output = new float[length];           //array of data with AGC applied
        int halfWindow = (int) Math.floor(window / 2.0);
        int start = halfWindow;                       //sample where sliding window begins
        int end = length - halfWindow - 1;            //sample where sliding window ends
        int nNonZero = 0;                             //count of non-zero samples in window
        double avg = 0;                               //average absolute value in window
        double abs = 0;                               //absolute value of the current sample
        double sum = 0;                               //sum of absolute values within window
        double defaultSFactor = Float.MAX_VALUE / 1000;

        if (window > length)
        {
            System.err.println("AGC: Invalid window length: " + window + " for trace length: " + length);
            return trace;
        }

        if (sFactor == 0.0)
        {
            System.err.println("AGC: Normalization Factor can't be zero!!! Being set to default value: " + defaultSFactor);
            sFactor = defaultSFactor;
        }

        //...Loop through samples within trace
        for (int iTrace = start; iTrace <= end; iTrace++)
        {
            nNonZero = 0;
            sum = 0;
            for (int iWindow = iTrace - halfWindow; iWindow <= iTrace + halfWindow; iWindow++)
            {
                abs = Math.abs(trace[iWindow]);
                if (agcDebug) {System.out.println("abs: " + abs);}
                if (abs > 0.0)
                {
                    nNonZero++;
                    sum += abs;
                }
            }
            avg = (nNonZero > 0) ?  sum / nNonZero : 0;
            scalars[iTrace] = (avg > 0) ? sFactor / avg : 1;   //if no non zero samples, do nothing to data
            if (agcDebug) {System.out.println("sum:" + sum + " avg: " + avg + " scalar: " + scalars[iTrace]);}
        }

        //...Set boundary condition of first and last samples getting the scalar from the nearest window
        for (int iTrace = 0; iTrace < start; iTrace++)
        {
            scalars[iTrace] = scalars[start];
        }
        for (int iTrace = end; iTrace < length; iTrace++)
        {
            scalars[iTrace] = scalars[end];
        }

        //...Apply AGC scalars
        double out = 0;
        for (int iTrace = 0; iTrace < length; iTrace++)
        {
            out = (scalars[iTrace] * (double) trace[iTrace]);
            output[iTrace] = (float) out;
            if (agcDebug)
            {
                System.out.println("data: " + trace[iTrace] + " scalar: " + scalars[iTrace] + " output: " + output[iTrace]);
            }
        }

        return output;
    }
    
    /**
     * AGC - Automatic Gain Control
     * Algorithm based on Yilmaz, Seismic Data Processing, 1987 page 61
     * Implementation mimics behavior of Paradigm's Focus AGC (handling of boundary conditions, ignoring zeroes, etc.)
     * <p/>
     * Summary:
     * Calculates average absolute amplitude for a given window.
     * The gain for that window will be = sFactor/avg_abs_amp.
     * The routine then moves the window up one sample and runs the calculation again.
     * <p/>
     * Once all of the gains are computed, data samples are multiplied by the gain value for each sample.
     * Samples at the beginning and end of the trace are given the gain from the nearest sample that
     * had a full window (i.e. - no partial windows are used).
     * <p/>
     * Error Handling:
     * If zero sFactor encountered, uses internally calculated sFactor (sets to reasonable value)
     * If window is greater than trace length, returns data unchanged.
     * 
     * Why Fast?
     * Recognizes that sums in middle of trace are the same except for adding one sample ahead of
     * trace and subtracting one sample behind trace. Significantly decreases looping.
     * 
     * NOTE:
     * window sizes should be odd (even sized windows are effectively rounded up to the nearest odd number)
     *
     * @param trace   Input trace samples (float array)
     * @param window  Length of sliding window (in samples)
     * @param sFactor Normalization amplitude (what you want the RMS amplitude of the output to be)
     * @return Output trace samples with AGC applied
     */
    public static float[] fastInstantaneousAGC(float[] trace, int window, double sFactor)
    {
        int length = trace.length;                    //length of input data array
        double[] scalars = new double[length];        //array off scalars to apply to data
        float[] output = new float[length];           //array of data with AGC applied
        int halfWindow = (int) Math.floor(window / 2.0);
        int start = halfWindow;                       //sample where sliding window begins
        int end = length - halfWindow - 1;            //sample where sliding window ends
        int nNonZero = 0;                             //count of non-zero samples in window
        double avg = 0;                               //average absolute value in window
        double abs = 0;                               //absolute value of the current sample
        double sum = 0;                               //sum of absolute values within window
        double defaultSFactor = Float.MAX_VALUE / 1000;

        if (window > length)
        {
            System.err.println("AGC: Invalid window length: " + window + " for trace length: " + length);
            return trace;
        }

        if (sFactor == 0.0)
        {
            System.err.println("AGC: Normalization Factor can't be zero!!! Being set to default value: " + defaultSFactor);
            sFactor = defaultSFactor;
        }

        //...Loop through beginning samples
        nNonZero = 0;
    	sum = 0;
        for (int iTrace=0; iTrace <= 2*halfWindow; iTrace++) {
        	abs = Math.abs(trace[iTrace]);
        	if (agcDebug) {System.out.println("abs: " + abs);}
        	if (abs > 0.0)
        	{
        		nNonZero++;
        		sum += abs;
        	}
        }
        avg = (nNonZero > 0) ?  sum / nNonZero : 0;
        scalars[start] = (avg > 0) ? sFactor / avg : 1;   //if no non zero samples, do nothing to data
        if (agcDebug) {System.out.println("sum: " + sum + " avg: " + avg + " scalar: " + scalars[start]);}
        
        //...Fill in beginning of scalar array
        for (int i=0;i<start;i++) {
        	scalars[i] = scalars[start];
        }
        
        //...Loop through middle traces
        for (int iTrace = start+1; iTrace <= end; iTrace++)
        {
        	//handle new sample
        	abs = Math.abs(trace[iTrace+halfWindow]);
        	if ( abs > 0.0 ) {
        		nNonZero++;
        		sum += abs;
        	}
        	if (agcDebug) {System.out.println("abs: " + abs + " sum: "+ sum + " nNonZero: "+ nNonZero);}
        	
        	//handle dropped sample
        	abs = Math.abs(trace[iTrace-halfWindow-1]);
        	if ( abs > 0.0 ) {
        		nNonZero--;
        		sum -= abs;
        	}
        	if (agcDebug) {System.out.println("abs: " + abs + " sum: "+ sum + " nNonZero: "+ nNonZero);}
        	
        	//compute scalar
        	avg = (nNonZero > 0) ?  sum / nNonZero : 0;
        	scalars[iTrace] = (avg > 0) ? sFactor / avg : 1;   //if no non zero samples, do nothing to data
        	if (agcDebug) {System.out.println("sum: " + sum + " avg: " + avg + " scalar: " + scalars[iTrace]);}
        }
    
    	//...Fill in end of scalar array
        for (int i=end;i<length;i++) {
        	scalars[i] = scalars[end];
        }

        //...Apply AGC scalars
        double out = 0;
        for (int iTrace = 0; iTrace < length; iTrace++)
        {
            out = (scalars[iTrace] * (double) trace[iTrace]);
            output[iTrace] = (float) out;
            if (agcDebug)
            {
                System.out.println("data: " + trace[iTrace] + " scalar: " + scalars[iTrace] + " output: " + output[iTrace]);
            }
        }

        return output;
    }

    public static float[] rmsAGC(float[] values, int halfWindow, double maxAmp)
    {
        int length = values.length;                    //length of input data array
        float[] output = new float[length];           //array of data with AGC applied                            //average absolute value in window
        double sumsq;                               //sumsq of absolute values within window
        double defaultSFactor = Float.MAX_VALUE / 1000;

        int window = 2 * halfWindow + 1;
        if (window > length)
        {
            System.err.println("AGC: Invalid window length: " + window + " for values length: " + length);
            return values;
        }

        if (maxAmp == 0.0)
        {
            System.err.println("AGC: Normalization Factor can't be zero!!! Being set to default value: " + defaultSFactor);
            maxAmp = defaultSFactor;
        }

        double[] vsq = new double[length];
        for (int i = 0; i < length; i++)
            vsq[i] = values[i] * values[i];

        double minAmpSq = 0.000001 * maxAmp * maxAmp;

        for (int i = 0; i < length; i++)
        {
            sumsq = 0;
            int start = Math.max(i - halfWindow, 0);
            int end = Math.min(i + halfWindow, length - 1);
            int windowSize = end - start + 1;
            for (int w = start; w <= end; w++)
            {
                sumsq += vsq[w];
            }
            double avgSq = sumsq / windowSize;
            output[i] = (avgSq <= minAmpSq) ? values[i] : (float) (values[i] / Math.sqrt(avgSq));
        }
        return output;
    }

    /** @param args  */
    public static void main(String[] args)
    {
        float[] data;
        float sFactor = 1;
        data = new float[]{0, 1, 0, -1, 0, 2, 0, -2, 0, 0, 1, 0, -1, 0, 2, 0, -2, 0, 0, 1, 0, -1, 0, 2, 0, -2, 0};
        
        int nAngles = 21;
        float[] sine = new float[nAngles];
        double dAngle = 360/nAngles;
        double angle = 0;
        for(int n = 0; n < nAngles; n++, angle += dAngle)
            sine[n] = (float) StsMath.sind(angle);

        StsToolkit.print(AGC.fastInstantaneousAGC(data, 6, 1), "fast");
        StsToolkit.print(AGC.instantaneousAGC(data, 6, 1), "slow");
/*
        testAGC(data, sFactor, "random data", 0, 4);
        timeAGC(data, sFactor, 4);
        /*
            StsMath.scale(data, 10);
            testAGC(data, sFactor, "random data scaled up by 10");

            data = new float[] { 0, 1, 0, -1, 0, 1, 0, -1, 0, 1, 0, -1, 0, 1, 0, -1, 0, 1, 0, -1, 0, 1, 0, -1, 0, 1, 0, -1, 0};
            testAGC(data, sFactor, "sinusoid data");
        */ /*
        int nAngles = 12;
        double[] sine = new double[nAngles];
        double dAngle = 360/nAngles;
        double angle = 0;
        for(int n = 0; n < nAngles; n++, angle += dAngle)
            sine[n] = StsMath.sind(angle);

        int nWaves = 10;
        int nn = 0;
        data = new float[nAngles*nWaves + 1];
        for(int n = 0; n < nWaves; n++)
            for(int a = 0; a < nAngles; a++, nn++)
                data[nn] = (float)sine[a];
        data[nn] = 0.0f;
        int halfWindowMin = nAngles/2;
        int halfWindowMax = 2*nAngles;
        testAGC(data, sFactor, nWaves + " sinusoid cycles in  " + nAngles + " divisions", halfWindowMin, halfWindowMax);
        timeAGC(data, sFactor, halfWindowMax);
        */
    }

    public static void testAGC(float[] data, float sFactor, String testName, int halfWindowMin, int halfWindowMax)
    {
        float[] out;
        StsToolkit.print(data, testName + "Input");
        int nSamples = data.length;
        System.out.println("\nComparison over range of window sizes\n");
        for (int halfWindow = halfWindowMin; halfWindow <= halfWindowMax; halfWindow++)
        {
            int window = 2 * halfWindow + 1;
            out = instantaneousAGC(data, 0, 0, nSamples, window, sFactor);
            rescale(out, sFactor);
            StsToolkit.print(out, "Focus instantaneousAGC window " + window);
            out = rmsAGC(data, 0, 0, nSamples, halfWindow, sFactor);
            rescale(out, sFactor);
            StsToolkit.print(out, "Focus rms AGC window " + window);
            out = StsSeismicFilter.applyAbsAGC(data, 0, 0, nSamples, halfWindow, 1.0f);
            rescale(out, sFactor);
            StsToolkit.print(out, "S2S abs AGC window " + window);
            out = StsSeismicFilter.applySuAGC(data, 0, 0, nSamples, halfWindow, 1.0f);
            rescale(out, sFactor);
            StsToolkit.print(out, "S2S rms AGC window " + window);
            System.out.println();
        }
    }

    private static float[] instantaneousAGC(float[] data, int window, double sFactor) {
    		if (fastAGC) return AGC.fastInstantaneousAGC(data, window, sFactor);
    		else return AGC.slowInstantaneousAGC(data, window, sFactor);
	}

	static private void rescale(float[] out, float maxAmp)
    {
        StsSeismicFilter.normalizeAmplitude(out);
        StsMath.scale(out, maxAmp);
    }

    public static void timeAGC(float[] data, float sFactor, int halfWindow)
    {
        int nSamples = data.length;
        int window = 2*halfWindow + 1;
        System.out.println("\nTiming comparison\n");
        StsTimer timer = new StsTimer();
        timer.start();
        for (int n = 0; n < 10000; n++)
            instantaneousAGC(data, 0, 0, nSamples, window, sFactor);
        double time = timer.stop();
        System.out.println("Focus Time per AGC (msec): " + time / 10);

        timer = new StsTimer();
        timer.start();
        for (int n = 0; n < 10000; n++)
            StsSeismicFilter.applyAbsAGC(data, 0, 0, nSamples, halfWindow, sFactor);
        time = timer.stop();
        System.out.println("S2S SU Time per AGC (msec): " + time/10);
    }
}
