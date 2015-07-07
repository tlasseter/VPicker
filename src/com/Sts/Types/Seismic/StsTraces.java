package com.Sts.Types.Seismic;

import com.Sts.UI.StsWiggleDisplayProperties;
import com.Sts.Utilities.StsMath;
import com.Sts.Utilities.Seismic.StsTraceUtilities;
import com.Sts.MVC.View3d.StsView2d;
import com.Sts.Types.*;

import javax.media.opengl.GL;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Nov 17, 2009
 * Time: 1:28:41 PM
 * To change this template use File | Settings | File Templates.
 */

/** This class handles traces to be isVisible in 2D: gathers, stacked sections, CVS, and VVS panels
 *  For each trace, equally-spaced subSamples may be interpolated.  If they are, isInterpolated is true.
 *  For a particular display situation, we may choose  to display interpolated data.  If we do, needsInterpolated is true.
 *  So the data at any time may or may not be interpolated.  If not, and we need interpolation, then we add the interpolated data.
 *  Note that we don't remove interpolated data once it has been computed.  We just decide whether we want to display it or not.
 *  So in the display loop, we may be skipping over interpolated data if we don't need to display it.
 *
 *  The data coming in is assumed to be normalized;  it's horizontal scale may be changed by various wiggleProperties.
 *  If we find that the horizontal scale has been changed, we compare this new scaling to the current scaling and multiply each value as needed.
 *  The new scaling is then saved as the current scaling.
 *
 *  This is an abstract class.  Concrete subclasses are for superGathers and stacked sections (StsTraceSet), and CVS and VVS panels (StsTracePanel).
 */

abstract public class StsTraces
{
    /** boundingBox for this data */
    StsRotatedGridBoundingBox boundingBox;
    /** number of traces */
    public int nTraces;
    /** number of original samples per trace */
    public int nSamples;
    /** number of equally-spaced subSamples */
    public int nInterpolatedSamples;
    /** minimum zValue for data */
    public float zMin;
   /** maximum zValue for data */
    public float zMax;
    /** zInc for data */
    public float zInc;
    /** subSampled zInc for gather data */
    public float zIncInterpolated;
    /** number of interpolated intervals from zInc to subSampled zInc */
    public int nInterpolationIntervals;
    /** zInc being isVisible (either zInc or zIncInterpolated */
    public float displayZInc;
    /** wiggle properties for this display */
    StsWiggleDisplayProperties wiggleProperties;
    /** Indicates values have not been interpolated yet */
    public boolean isInterpolated = false;
    /** Indicates trace should be drawn without interpolated points */
    boolean needsInterpolated = false;
    /** current scale value being used */
    double currentScale = 1.0;
    /** sample display range min.  If interpolated, this is the index into the interpolated array. */
    int sampleDisplayMin;
    /** sample display range max.  If interpolated, this is the index into the interpolated array. */
    int sampleDisplayMax;
    /** sample increment */
    int sampleDisplayInc;
    /** number of samples to display. Takes into consideration whether values are interpolated and interpolated values need to be isVisible. */
    int nDisplaySamples;
    /** z value of first sample */
    float sampleZMin;

    static public final boolean debug = false;

    abstract public void interpolateData();
    abstract public void scaleData(float scaleFactor);


    public StsTraces()
    {
    }

    public StsTraces(StsRotatedGridBoundingBox boundingBox)
    {
        this.boundingBox = boundingBox;
        this.zMin = boundingBox.zMin;
        this.zMax = boundingBox.zMax;
        this.zInc = boundingBox.zInc;
        this.nSamples = boundingBox.nSlices;
    }

    public StsTraces(StsRotatedGridBoundingBox boundingBox, int nApproxIntervals)
    {
        this(boundingBox);
        initializeInterpolation(nApproxIntervals);
    }

    private void initializeInterpolation(int nApproxIntervals)
    {
        this.nInterpolationIntervals = StsTraceUtilities.computeInterpolationInterval(zInc, nApproxIntervals);
        zIncInterpolated = zInc/nInterpolationIntervals;
        nInterpolatedSamples = (nSamples - 1)*nInterpolationIntervals + 1;
    }

    public void initializeDraw(StsView2d view2d, StsWiggleDisplayProperties wiggleProperties)
    {
        initializeDraw(view2d, 1.0, wiggleProperties);   
    }

    public void initializeDraw(StsView2d view2d, double horizScaleFactor, StsWiggleDisplayProperties wiggleProperties)
    {
        this.wiggleProperties = wiggleProperties;
		// TODO FIX: should only have to do this when view changed or moved
        view2d.computeVerticalPixelScaling();
		float pixelsPerYunit = Math.abs(view2d.pixelsPerYunit);
        float[][] axisRanges = view2d.axisRanges;
        float[] displayZRange = StsTraceUtilities.getDisplayZRange(axisRanges, zInc);
        int[] sampleDisplayRange = StsTraceUtilities.getDisplayRange(boundingBox, displayZRange);
        int pixelsPerInc = StsMath.ceiling(pixelsPerYunit * zInc / 2);
        needsInterpolated = wiggleProperties.getWiggleSmoothCurve();
        if (pixelsPerInc <= 2) needsInterpolated = false;
        if(needsInterpolated && !isInterpolated)
        {
            interpolateData();
            isInterpolated = true;
        }
        float overlapPercent = wiggleProperties.getWiggleOverlapPercent();
        double scale = horizScaleFactor*(100 + overlapPercent)/100;
        //  double scale = horizScaleFactor*(overlapPercent/100;
        boolean reversePolarity = wiggleProperties.getWiggleReversePolarity();
        if (reversePolarity) horizScaleFactor = -horizScaleFactor;
        if(scale != currentScale)
        {
            scaleData((float)(scale/currentScale));
            currentScale = scale;
        }

        sampleDisplayMin = sampleDisplayRange[0];
        sampleDisplayMax = sampleDisplayRange[1];
        sampleZMin = zMin + sampleDisplayMin*zInc;
        sampleDisplayInc = 1;
        displayZInc = zInc;
        nDisplaySamples = sampleDisplayMax - sampleDisplayMin + 1;
        if(isInterpolated)
        {
            if(needsInterpolated)
            {
                sampleDisplayInc = 1;
                displayZInc = zIncInterpolated;
                nDisplaySamples = (sampleDisplayMax - sampleDisplayMin)*nInterpolationIntervals + 1;
            }
            else
            {
                sampleDisplayInc = nInterpolationIntervals;
                displayZInc = zInc;
            }
            sampleDisplayMin *= nInterpolationIntervals;
            sampleDisplayMax *= nInterpolationIntervals;
        }

    }

    public void displayInterpolatedPoints(GL gl, float[] trace, float y, StsWiggleDisplayProperties wiggleProperties)
    {
        if(wiggleProperties.hasFill())
        {
            StsTraceUtilities.drawFilledWiggleTraces(gl, trace, sampleDisplayMin, nDisplaySamples, y, sampleZMin, displayZInc, wiggleProperties, sampleDisplayInc);
            StsTraceUtilities.drawFilledWiggleTracesLine(gl, trace, sampleDisplayMin, nDisplaySamples, y, sampleZMin, displayZInc, wiggleProperties, sampleDisplayInc);
        }
        if (wiggleProperties.getWiggleDrawLine())
            StsTraceUtilities.drawWiggleTraces(gl, trace, sampleDisplayMin, nDisplaySamples, y, sampleZMin, displayZInc, wiggleProperties, sampleDisplayInc);
        if (wiggleProperties.getWiggleDrawPoints())
            StsTraceUtilities.drawWigglePoints(gl, trace, sampleDisplayMin, nDisplaySamples, y, sampleZMin, displayZInc, sampleDisplayInc); // uses tracePoints and tracePointTypes in StsTraceUtilities
    }

    public void displayInterpolatedPoints(GL gl, float[] trace, float y, double[] muteRange, StsWiggleDisplayProperties wiggleProperties)
    {
        if(wiggleProperties.hasFill())
        {
            StsTraceUtilities.drawFilledWiggleTraces(gl, trace, sampleDisplayMin, nDisplaySamples, y, sampleZMin, displayZInc, wiggleProperties, muteRange, sampleDisplayInc);
            StsTraceUtilities.drawFilledWiggleTracesLine(gl, trace, sampleDisplayMin, nDisplaySamples, y, sampleZMin, displayZInc, wiggleProperties, muteRange, sampleDisplayInc);
        }
        if (wiggleProperties.getWiggleDrawLine())
            StsTraceUtilities.drawWiggleTraces(gl, trace, sampleDisplayMin, nDisplaySamples, y, sampleZMin, displayZInc, wiggleProperties, sampleDisplayInc);
        if (wiggleProperties.getWiggleDrawPoints())
            StsTraceUtilities.drawWigglePoints(gl, trace, sampleDisplayMin, nDisplaySamples, y, sampleZMin, displayZInc, sampleDisplayInc); // uses tracePoints and tracePointTypes in StsTraceUtilities
    }
}
