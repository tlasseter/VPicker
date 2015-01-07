package com.Sts.Types.Seismic;

import com.Sts.UI.StsWiggleDisplayProperties;
import com.Sts.Types.*;

import javax.media.opengl.GL;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Nov 18, 2009
 * Time: 2:40:42 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsTraceSet extends StsTraces
{
    /** data for traces */
    public StsTrace[] traces;

    public StsTraceSet()
    {
    }

    public StsTraceSet(StsRotatedGridBoundingBox boundingBox)
    {
        super(boundingBox);
    }

    public StsTraceSet(StsRotatedGridBoundingBox boundingBox, int nApproxInterpIntervals)
    {
        super(boundingBox, nApproxInterpIntervals);
    }

    public StsTraceSet(float[][] data, StsRotatedGridBoundingBox boundingBox, int nInterpolationIntervals, StsWiggleDisplayProperties wiggleProperties)
    {
        super(boundingBox, nInterpolationIntervals);
        nTraces = data.length;
        nSamples = data[0].length;

        traces = new StsTrace[nTraces];
        for(int n = 0; n < nTraces; n++)
            traces[n] = new StsTrace(data[n]);
    }

    public void interpolateData()
    {
        if(isInterpolated) return;
        for(int n = 0; n < nTraces; n++)
            traces[n].interpolateData(nInterpolationIntervals, zInc);
        isInterpolated = true;
    }

    public void computePointTypes()
    {
        for(int n = 0; n < nTraces; n++)
            traces[n].computePointTypes();
    }

    public void scaleData(float scaleMultiplier)
    {
        for(int n = 0; n < nTraces; n++)
            traces[n].scaleData(scaleMultiplier);
    }

    public void displayInterpolatedPoints(GL gl, int traceIndex, float y)
    {
        displayInterpolatedPoints(gl, traces[traceIndex].data, y, wiggleProperties);
    }

    public void displayInterpolatedPoints(GL gl, int traceIndex, float y, double[] muteRange)
    {
        displayInterpolatedPoints(gl, traces[traceIndex].data, y, muteRange, wiggleProperties);
    }

    public void displayInterpolatedPoints(GL gl, StsTrace trace, float y, double[] muteRange)
    {
        displayInterpolatedPoints(gl, trace.data, y, muteRange, wiggleProperties);
    }

    public synchronized boolean computeInstantaneousAmpAndPhaseGatherData()
    {
        interpolateData();
        computePointTypes();
        for (int n = 0; n < nTraces; n++)
//            traces[n].computeInstantaneousAmpAndPhaseGatherData(nInterpolationIntervals, nInterpolatedSamples, zInc);
        traces[n].computeInstantaneousAmpAndPhaseGatherData();       // jbw merge problem - i dunno

        return true;
    }

    public double computeRmsAmplitude()
    {
        double rmsAmplitude = 0.0;
        for(StsTrace trace : traces)
            rmsAmplitude = Math.max(rmsAmplitude, trace.computeRmsAmplitude());
        return rmsAmplitude;
    }

    public void normalizeAmplitude(double amplitude)
    {
        for(StsTrace trace : traces)
            trace.normalizeAmplitude(amplitude);
    }

    public StsTrace getFirstTrace() { return traces[0]; }
    public StsTrace getLastTrace() { return traces[nTraces-1]; }
}
