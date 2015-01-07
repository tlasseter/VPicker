package com.Sts.DBTypes;

/**
 * Title:        Sts Well Viewer Description:  Well Model-Viewer Copyright:    Copyright (c) 2001
 * Company:      4D Systems LLC
 * @author T.J.Lasseter
 * @version 1.0
 */

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.MVC.ViewWell.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.Utilities.*;

public class StsLogCurveType extends StsMainObject implements StsTreeObjectI
{
    protected StsObjectRefList logCurves;
    protected float curveMin = StsParameters.largeFloat;
    protected float curveMax = -StsParameters.largeFloat;
    protected float displayCurveMin;
    protected float displayCurveMax;
    protected int scaleType = LIN;
    protected StsLogCurveType aliasToType = this;
    /** index of color in standard spectrum used in drawing log curve line */
    protected int colorIndex;

    // transient protected String aliasToString = STRING_NONE;
    transient double[] scale;
    transient StsLogCurveTypeClass logCurveTypeClass = null;

    final static public int LIN = 0;
    final static public int LOG = 1;

    final static public String[] scaleTypeStrings = new String[] { "LIN", "LOG" };

    static public String STRING_NONE = "None";

    static public StsFieldBean[] displayFields = null;
    static StsObjectPanel objectPanel = null;

    public StsLogCurveType()
    {
    }

    public StsLogCurveType(boolean persistent)
    {
        super(persistent);
    }

    public StsLogCurveType(StsLogCurve logCurve)
    {
		super(false);
        String name = logCurve.getName();
        setName(name);
        initialize();
        colorIndex = logCurveTypeClass.getNextColorIndex(name);
		addToModel();
		adjustRange(logCurve);
		logCurves = StsObjectRefList.constructor(2, 2, "logCurves", this);
        initializeScale();
    }

    public StsLogCurveType(String name, float min, float max)
    {
		super(false);
        setName(name);
        curveMin = min;
        curveMax = max;
        initialize();
        colorIndex = logCurveTypeClass.getNextColorIndex(name);
		logCurves = StsObjectRefList.constructor(2, 2, "logCurves", this);
        initializeScale();
        addToModel();
    }

    private void initialize()
    {
        logCurveTypeClass = (StsLogCurveTypeClass)getStsClass();
    }

    public boolean initialize(StsModel model)
    {
        initialize();
    /*
        if(aliasToType == null)
            aliasToString = STRING_NONE;
        else
            aliasToString = aliasToType.getName();
    */
        reinitializeScale();
        return true;
    }

    public boolean hasAlias()
    {
        return aliasToType != this;
    }
/*
    public String getAliasToString()
    {
        return aliasToType.getName();
    }

    public void setAliasToString(String aliasToString)
    {
        if(this.aliasToString.equals(aliasToString))
            this.aliasToString = STRING_NONE;
        else
            this.aliasToString = aliasToString;

        if(aliasToString == STRING_NONE)
            aliasToType = null;
        else
        {
            aliasToType = (StsLogCurveType)currentModel.getObjectWithName(StsLogCurveType.class, aliasToString);
            if(aliasToType != null)
            {
                dbFieldChanged("aliasToType", aliasToType);
                colorIndex = aliasToType.getColorIndex();
                setDisplayCurveMin(aliasToType.getDisplayCurveMin());
                setDisplayCurveMax(aliasToType.getDisplayCurveMax());
            }
        }
        viewChanged();
    }
*/
    static public String[] getLogCurveTypeStrings(StsModel model)
    {
        return model.getInstanceNames(StsLogCurveType.class);
    }

    private void viewChanged()
    {
        currentModel.repaintViews(StsView3d.class);
        currentModel.repaintViews(StsLogCurvesView.class);
    }

    public StsColor getStsColor() { return StsLogCurveTypeClass.logCurveTypeColors[colorIndex]; }

    public void setStsColor(StsColor stsColor)
    {
        int index = StsColor.getColorIndex(stsColor, StsLogCurveTypeClass.logCurveTypeColors);
        if(index == colorIndex) return;
        colorIndex = index;
        dbFieldChanged("colorIndex", index);
        viewChanged();
    }

    public int getColorIndex() { return colorIndex; }

    public float getCurveMin() { return curveMin; }
    public float getCurveMax() { return curveMax; }

    public float getDisplayCurveMin() { return displayCurveMin; }
    public void setDisplayCurveMin(float min)
    {
        if(displayCurveMin == min) return;
        displayCurveMin = min;
        dbFieldChanged("displayCurveMin", min);
        scale = StsMath.niceScale(displayCurveMin, displayCurveMax, 10, isLinear());
        viewChanged();
    }

    public float getDisplayCurveMax() { return displayCurveMax; }
    public void setDisplayCurveMax(float max)
    {
        if(displayCurveMax == max) return;
        displayCurveMax = max;
        dbFieldChanged("displayCurveMax", max);
        scale = StsMath.niceScale(displayCurveMin, displayCurveMax, 10, isLinear());
        viewChanged();
    }

    public void setDisplayScale(float min, float max)
    {
        displayCurveMin = min;
        displayCurveMax = max;
    }

    public boolean isLinear() { return scaleType == LIN; }

    public double[] getScale() { return scale; }

    public double getScaleMin() { return scale[0]; }
    public double getScaleMax() { return scale[1]; }

    public String getScaleTypeString() { return scaleTypeStrings[scaleType]; }
    public void setScaleTypeString(String scaleTypeString)
    {
        for(int n = 0; n < 2; n++)
        {
            if(scaleTypeString == scaleTypeStrings[n])
            {
                scaleType = n;
                dbFieldChanged("scaleType", scaleType);
                scale = StsMath.niceScale(displayCurveMin, displayCurveMax, 10, isLinear());
                viewChanged();
            }
        }
    }

    public int getScaleType() { return scaleType; }

    protected void adjustRange(StsLogCurve logCurve)
    {
        float logCurveMin = logCurve.getMinValue();
        float logCurveMax = logCurve.getMaxValue();
        adjustRange(logCurveMin, logCurveMax);
    }

    protected void adjustRange(float logCurveMin, float logCurveMax)
    {
        if(logCurveMin < curveMin)
        {
            curveMin = logCurveMin;
            dbFieldChanged("curveMin", curveMin);
        }
        if(logCurveMax > curveMax)
        {
            curveMax = logCurveMax;
            dbFieldChanged("curveMax", curveMax);
        }
		initializeScale();
    }

	protected void initializeScale()
	{
        scaleType = computeScaleType(curveMin, curveMax);
        scale = StsMath.niceScale(curveMin, curveMax, 10, isLinear());
        displayCurveMin = (float)scale[0];
        displayCurveMax = (float)scale[1];
    }

	protected void reinitializeScale()
	{
        scale = StsMath.niceScale(displayCurveMin, displayCurveMax, 10, isLinear());
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

    static public int computeScaleType(StsLogCurve logCurve)
    {
        return computeScaleType(logCurve.getMinValue(), logCurve.getMaxValue());
    }

    static public int computeScaleType(float curveMin, float curveMax)
    {
        if (curveMin > 0.0f && curveMax / curveMin > 100.0f)
            return LOG;
        else
            return LIN;
    }

    public void addLogCurve(StsLogCurve logCurve)
    {
        if (logCurve == null) return;
		if (logCurves == null)
	    {
			logCurves = StsObjectRefList.constructor(2, 2, "logCurves", this);
		}
        logCurves.add(logCurve);
		adjustRange(logCurve);
    }

    static public Object[] getLogCurveTypes()
    {
        return currentModel.getObjectList(StsLogCurveType.class);
    }

    public StsFieldBean[] getDisplayFields()
    {
        if(displayFields == null)
        {
            displayFields = new StsFieldBean[]
            {
                new StsStringFieldBean(StsLogCurveType.class, "name", false, "Name:"),
                new StsComboBoxFieldBean(StsLogCurveType.class, "aliasToType", "Alias to:", "logCurveTypes"),
                new StsFloatFieldBean(StsLogCurveType.class, "displayCurveMin", true, "Display min:"),
                new StsFloatFieldBean(StsLogCurveType.class, "displayCurveMax", true, "Display max:"),
                new StsFloatFieldBean(StsLogCurveType.class, "curveMin", false, "Curve min:"),
                new StsFloatFieldBean(StsLogCurveType.class, "curveMax", false, "Curve max:"),
                new StsColorComboBoxFieldBean(StsLogCurveType.class, "stsColor", "Color:", StsLogCurveTypeClass.logCurveTypeColors),
                new StsComboBoxFieldBean(StsLogCurveType.class, "scaleTypeString", "Curve Type:", scaleTypeStrings)
            };
        }
        return displayFields;
    }

    public StsFieldBean[] getPropertyFields() { return null; }

    public StsFieldBean[] getDefaultFields() { return null; }

    public Object[] getChildren()
    {
        return new Object[0];
    }

    public boolean anyDependencies()
    {
        return false;
    }

    public StsObjectPanel getObjectPanel()
    {
        if (objectPanel == null)
        {
            objectPanel = StsObjectPanel.constructor(this, true);
        }
        return objectPanel;
    }

    public void treeObjectSelected()
    {
        getStsClass().selected(this);
    }

}
