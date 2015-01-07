package com.Sts.UI;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.Utilities.*;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */
public abstract class StsSeismicPanelProperties extends StsPanelProperties implements StsSerializable
{
    // CVS & Semblance
    public float minVelocityStep = 100.0f;
    public float velocityStep = 0.1f;
    public float velocityMin = StsPreStackVelocityModel.DefaultVelMin;
    public float velocityMax = StsPreStackVelocityModel.DefaultVelMax;

    public int wiggleOverlapPercent;

    // CVS, Semblance, Backbone
    public boolean showLabels = true;

    /** checked by all appropriate views on a viewObjectChanged(seismicPanelProperties) to see if inchesPerSecond or tracesPerInch has changed */
    private transient boolean rangeChanged = false;

    public StsSeismicPanelProperties()
	{
	}

	public StsSeismicPanelProperties(String panelTitle, String fieldName)
	{
        super(panelTitle, fieldName);
        initializeVelMaxMin();
	}

	public StsSeismicPanelProperties(Object parentObject, Object defaultProperties, String panelTitle, String fieldName)
	{
        super(parentObject, panelTitle, fieldName);
        initializeDefaultProperties(defaultProperties);
        initializeVelMaxMin();
	}
	
    public void initializeVelMaxMin()
    {
        StsPreStackLineSet lineSet = StsPreStackLineSetClass.currentProjectPreStackLineSet;
        if (lineSet != null && lineSet.velocityModel != null)
        {
            velocityMin = lineSet.velocityModel.dataMin;
            velocityMax = lineSet.velocityModel.dataMax;
        }
    }

    public void initializeDefaultProperties(Object defaultProperties)
    {
        if(defaultProperties != null)
            StsToolkit.copyObjectNonTransientFields(defaultProperties, this);
    }
    // Getters
    public boolean getShowLabels() { return showLabels; }
/*
    public int getCorridorPercentage() { return corridorPercentage; }
    public boolean getShowAvgProfile() { return showAvgProfile; }
    public boolean getShowNeighborProfiles() { return showNeighborProfiles; }
    public boolean getShowLimitProfiles() { return showLimitProfiles; }
    public boolean getShowIntervalVelocities() { return showIntervalVelocities; }
*/
    public float getVelocityMin() { return velocityMin * 1000; }
    public float getVelocityMax() { return velocityMax * 1000; }
    public float getVelocityStep() { return velocityStep * 1000; }
    public int getWiggleOverlapPercent() { return wiggleOverlapPercent; }
 //   public float getTracesPerInch() { return tracesPerInch; }
 //   public float getInchesPerSecond() { return inchesPerSecond; }

    // Setters
    public void setVelocityMin(float min)
    { 
    	velocityMin = min/1000; rangeChanged = true;
    	resetVelocityColorscale();
    }
    public void setVelocityMax(float max) 
    { 
    	velocityMax = max/1000; rangeChanged = true;
    	resetVelocityColorscale();
    }
    
    private void resetVelocityColorscale() 
    {
        StsPreStackLineSet lineSet = StsPreStackLineSetClass.currentProjectPreStackLineSet;
        if (lineSet != null) { 
            if (lineSet.velocityModel != null) {
                lineSet.velocityModel.setVelocityMaxMin(velocityMax, velocityMin);
            }
        }
    }
    
    public void setVelocityStep(float step)
    {
        if(step <= 0.0) return;
        velocityStep = step / 1000;
        rangeChanged = true;
    }
    
    public void setShowLabels(boolean showLabels) 
    { 
        this.showLabels = showLabels; 
    }
    
    public void setWiggleOverlapPercent(int percent) 
    { 
        wiggleOverlapPercent = percent; 
    }
/*
    public void setCorridorPercentage(int percentage) { corridorPercentage = percentage; }
    public void setShowAvgProfile(boolean show) { this.showAvgProfile = show; }
    public void setShowNeighborProfiles(boolean show) { this.showNeighborProfiles = show; }
    public void setShowLimitProfiles(boolean show) { this.showLimitProfiles = show; }
    public void setShowIntervalVelocities(boolean show) { this.showIntervalVelocities = show; }
    public void setTracesPerInch(Float tracesPerInch) { setTracesPerInch(tracesPerInch.floatValue()); }
    public void setTracesPerInch(float tracesPerInch)
    {
        if(this.tracesPerInch == tracesPerInch) return;
        this.tracesPerInch = tracesPerInch;
        rangeChanged = true;
//        currentModel.viewObjectChanged(this);
//        currentModel.viewObjectRepaint(this);
    }
    public void setInchesPerSecond(Float inchesPerSecond) { setInchesPerSecond(inchesPerSecond.floatValue()); }
    public void setInchesPerSecond(float inchesPerSecond)
    {
        if(this.inchesPerSecond == inchesPerSecond)return;
        this.inchesPerSecond = inchesPerSecond;
        rangeChanged = true;
//        currentModel.viewObjectChanged(this);
//        currentModel.viewObjectRepaint(this);
    }
*/
	public void dialogSelectionType(int selectionType)
	{
        super.dialogSelectionType(selectionType);
        rangeChanged = false;
    }

    public boolean isRangeChanged()
    {
        return rangeChanged;
    }

    public void setRangeChanged(boolean rangeChanged)
    {
        this.rangeChanged = rangeChanged;
    }
}
