package com.Sts.DBTypes;

/**
 * Title:        Sts Well Viewer Description:  Well Model-Viewer Copyright:    Copyright (c) 2001
 * Company:      4D Systems LLC
 * @author T.J.Lasseter
 * @version 1.0
 */

import com.Sts.MVC.*;
import com.Sts.Utilities.*;

public class StsTimeCurveType extends StsObject
{
    protected String name;
    protected StsObjectRefList timeCurves;
    protected float curveMin = StsParameters.largeFloat;
    protected float curveMax = -StsParameters.largeFloat;
    transient float displayCurveMin;
    transient float displayCurveMax;
    transient int displayScaleType = LIN;
    transient double[] scale;
    protected StsColor stsColor = null;
    final static public int LIN = 0;
    final static public int LOG = 1;

    final static public String[] scaleTypeLabels = new String[] { "LIN", "LOG" };

    static private String[] curveNames = { "PRES", "TEMP", "SATUR", "WCUT" };
    static private StsColor[] curveColors = { StsColor.YELLOW, StsColor.GREEN, StsColor.CYAN, StsColor.BLUE};
    static StsSpectrum colorSpectrum;

    public StsTimeCurveType()
    {
    }

    public StsTimeCurveType(StsTimeCurve timeCurve, StsColor color)
    {
		super(false);
        this.name = timeCurve.getName();

        if (color != null)
            stsColor = color;
        else
            stsColor = incrementCurrentColor();
		adjustMinMax(timeCurve);
        addToModel();
        timeCurves = StsObjectRefList.constructor(2, 2, "timeCurves", this);
    }

    public boolean initialize(StsModel model)
    {
        if(timeCurves == null)
            timeCurves = StsObjectRefList.constructor(2, 2, "timeCurves", this);

		setCurveDisplay();
        return true;
    }

    private StsColor incrementCurrentColor()
    {
        if(colorSpectrum == null)
        {
            colorSpectrum = new StsSpectrum("Basic");
            colorSpectrum.setBasic32Colors();
			colorSpectrum.addToModel();
		}
        return colorSpectrum.incrementCurrentColor();
    }

    public String getName() { return name; }

    public StsColor getStsColor() { return stsColor; }

    public float getDisplayCurveMin() { return displayCurveMin; }
    public float getDisplayCurveMax() { return displayCurveMax; }

    public boolean isLinear() { return displayScaleType == LIN; }

    public double[] getScale() { return scale; }

    public double getScaleMin() { return scale[0]; }
    public double getScaleMax() { return scale[1]; }

    public String getScaleTypeLabel() { return scaleTypeLabels[displayScaleType]; }

    public int getScaleType() { return displayScaleType; }

    private void adjustMinMax(StsTimeCurve timeCurve)
    {
        setCurveMin(timeCurve.getCurveMin());
        setCurveMax(timeCurve.getCurveMax());
		setCurveDisplay();
    }

	private void setCurveDisplay()
	{
        displayScaleType = computeScaleType(curveMin, curveMax);
        scale = StsMath.niceScale(curveMin, curveMax, 10, isLinear());
        displayCurveMin = (float)scale[0];
        displayCurveMax = (float)scale[1];
    }

	private void setCurveMin(float value)
	{
		curveMin = Math.min(curveMin, value);
		dbFieldChanged("curveMin", value);
	}

	private void setCurveMax(float value)
	{
		curveMax = Math.max(curveMax, value);
		dbFieldChanged("curveMax", value);
	}

    static public int computeScaleType(StsTimeCurve timeCurve)
    {
        return computeScaleType(timeCurve.getMinValue(), timeCurve.getMaxValue());
    }

    static public int computeScaleType(float curveMin, float curveMax)
    {
        if (curveMin > 0.0f && curveMax / curveMin > 100.0f)
            return LOG;
        else
            return LIN;
    }

    public void addLogCurve(StsTimeCurve timeCurve)
    {
        if (timeCurve == null) return;
        if(timeCurves == null)
            timeCurves = StsObjectRefList.constructor(2, 2, "timeCurves", this);
        timeCurves.add(timeCurve);
		adjustMinMax(timeCurve);
    }
}
