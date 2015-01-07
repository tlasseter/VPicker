package com.Sts.UI;

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

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Types.PreStack.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

public class StsVVSProperties extends StsPanelProperties
{
    public float topPercent = 3.0f;
    public float btmPercent = 4.0f;

    protected int numberPanels = 5;
	protected boolean ignoreSuperGather = false;
    protected int tracesPerPanel = 5;
    protected int wiggleOverlapPercent = 30;
    protected float inchesPerSecond = 3.f;

    public boolean showLabels = true;
    public boolean showAvgProfile = true;
    public boolean showNeighborProfiles = true;
    public boolean showLimitProfiles = true;
    public boolean showIntervalVelocities = true;
    public int ilineDisplayInc = 10;
    public int xlineDisplayInc = 10;
    public int corridorPercentage = 100;

//    protected float tracesPerInch = 15.f;
//    protected transient float currentTracesPerInch = tracesPerInch;
    protected transient float currentInchesPerSecond = inchesPerSecond;

	transient StsGroupBox groupBox = null;
    transient StsBooleanFieldBean isCVSBean;
    transient StsFloatFieldBean topPercentBean;
    transient StsFloatFieldBean btmPercentBean;

	transient StsFloatFieldBean velocityStepBean;
    transient StsFloatFieldBean velocityMinBean;
    transient StsFloatFieldBean velocityMaxBean;
    transient StsIntFieldBean numberPanelsBean;
    transient StsIntFieldBean tracesPerPanelBean;
	transient StsBooleanFieldBean ignoreSuperGatherBean;
//    transient StsFloatFieldBean tpiBean;
    transient StsFloatFieldBean ipsBean;
    transient StsIntFieldBean wiggleOverlapPercentBean;

	transient public StsModel model;
	transient public boolean recompute = true;

	static private final String title = "VVS Display Properties";

	public StsVVSProperties()
	{
	}

	public StsVVSProperties(StsModel model, StsClass stsClass, String fieldName)
	{
		super(title, fieldName);
		this.model = model;
//        ilineDisplayInc = StsPreStackLineSet.defaultAnalysisRowInc;
//        xlineDisplayInc = StsPreStackLineSet.defaultAnalysisColInc;
//        corridorPercentage = StsPreStackLineSet.defaultCorridorPercentage;
	}

	public StsVVSProperties(StsObject parentObject, StsVVSProperties defaultProperties, String fieldName)
	{
        super(parentObject, "VVS Properties", fieldName);
        initializeDefaultProperties(defaultProperties);
    }
    public void initializeDefaultProperties(Object defaultProperties)
    {
        if(defaultProperties != null)
            StsToolkit.copyObjectNonTransientFields(defaultProperties, this);
    }

    public void initializeVelocity(StsObject parentObject)
    {
        if(!ignoreSuperGather)
           tracesPerPanel = getNGathersInSuperGather();
    }

	public void initializeBeans()
	{
        tracesPerPanel = getNGathersInSuperGather();
        wiggleOverlapPercentBean = new StsIntFieldBean(this, "wiggleOverlapPercent", 0, 1000, "Gain (Overlap Percent):", true);
        wiggleOverlapPercentBean.setStep(25);
		propertyBeans = new StsFieldBean[]
		{
//            tpiBean = new StsFloatFieldBean(this, "tracesPerInch", 1.f, 250.f, "Traces per Inch:"),
            topPercentBean = new StsFloatFieldBean(this, "topPercent", 0, 100, "Top Percentage:", true),
            btmPercentBean = new StsFloatFieldBean(this, "btmPercent", 0, 100, "Bottom Percentage:", true),
            ignoreSuperGatherBean = new StsBooleanFieldBean(this, "ignoreSuperGather", "Ignore Super Gather:"),
            tracesPerPanelBean = new StsIntFieldBean(this, "tracesPerPanel", 0, 10, "Number of Traces/Panel:", true),
            numberPanelsBean = new StsIntFieldBean(this, "numberPanels", 0, 20, "Number of Panels:", true),

            wiggleOverlapPercentBean,
            ipsBean = new StsFloatFieldBean(this, "inchesPerSecond", .1f, 50.0f, "Inches per Second:"),
            new StsIntFieldBean(this, "ilineDisplayInc", 0, 100, "Inline Profile Display Distance:", true),
            new StsIntFieldBean(this, "xlineDisplayInc", 0, 100, "Crossline Profile Display Distance:", true),
            new StsBooleanFieldBean(this, "showLabels", "Show Point Labels"),
            new StsBooleanFieldBean(this, "showAvgProfile", "Show Avg Profile"),
            new StsBooleanFieldBean(this, "showNeighborProfiles", "Show Neighbor Profiles"),
            new StsBooleanFieldBean(this, "showLimitProfiles", "Show Limit Profiles"),
            new StsIntFieldBean(this, "corridorPercentage", 10, 100, "Corridor Percentage:", true),
            new StsBooleanFieldBean(this, "showIntervalVelocities", "Show Interval Velocities"),
            };
        velocityStepBean.fixStep(100);
        topPercentBean.fixStep(.1);
        btmPercentBean.fixStep(.1);

        if(!ignoreSuperGather)
            tracesPerPanelBean.setEditable(false);
	}

    private int getNGathersInSuperGather()
    {
        StsSuperGather superGather = ((StsPreStackLineSet)parentObject).getSuperGather();
        if(superGather == null) return 0;
        else return superGather.getNumberGathers();
    }
    public int getNumberPanels() { return numberPanels; }
	public boolean getIgnoreSuperGather() { return ignoreSuperGather; }
    public int getTracesPerPanel() { return tracesPerPanel; }

    public void setNumberPanels(int nPanels)
    {
        if(nPanels <= 0)
            return;
        numberPanels = nPanels;
    }
    public void setIgnoreSuperGather(boolean use)
    {
        ignoreSuperGather = use;
        if(ignoreSuperGather)
            tracesPerPanelBean.setEditable(true);
        else
        {
            tracesPerPanelBean.setEditable(false);
            tracesPerPanel = getNGathersInSuperGather();
        }
    }

    public void setWiggleOverlapPercent(int percent)
    {
        wiggleOverlapPercent = percent;
	}
    public int getWiggleOverlapPercent()
    {
        return wiggleOverlapPercent;
	}
    public float getTopPercent()
    {
        return topPercent;
    }
    public void setTopPercent(float percent)
    {
        topPercent = percent;
    }
    public float getBtmPercent()
    {
        return btmPercent;
    }
    public void setBtmPercent(float percent)
    {
        btmPercent = percent;
    }

    public void setTracesPerPanel(int ntraces) { tracesPerPanel = ntraces; }
/*
    public void setTracesPerInch(Float tracesPerInch)
    {
        setTracesPerInch(tracesPerInch.floatValue());
    }

    public void setTracesPerInch(float tracesPerInch)
    {
        setTracesPerInch(tracesPerInch, true);
    }

    public void setTracesPerInch(float tracesPerInch, boolean setval)
    {
        this.currentTracesPerInch = tracesPerInch;
        if(tpiBean != null) tpiBean.setValue(tracesPerInch);
        if(!setval) return;
        if(this.tracesPerInch == tracesPerInch)return;
        this.tracesPerInch = tracesPerInch;
        currentModel.viewObjectChanged(this);
        currentModel.viewObjectRepaint(this);
    }
*/
    public void setInchesPerSecond(Float inchesPerSecond)
    {
        setInchesPerSecond(inchesPerSecond.floatValue());
    }

    public void setInchesPerSecond(float inchesPerSecond)
    {
        setInchesPerSecond(inchesPerSecond, true);
    }

    public void setInchesPerSecond(float inchesPerSecond, boolean setval)
    {
        this.currentInchesPerSecond = inchesPerSecond;
        if(ipsBean != null) ipsBean.setValue(inchesPerSecond);
        if(!setval) return;
        if(this.inchesPerSecond == inchesPerSecond)return;
        this.inchesPerSecond = inchesPerSecond;
        currentModel.viewObjectChanged(this, this);
        currentModel.viewObjectRepaint(this, this);
    }
/*
    public float getTracesPerInch()
    {
        return tracesPerInch;
    }
*/
    public float getInchesPerSecond()
    {
        return inchesPerSecond;
    }
    public boolean getShowLabels() { return showLabels; }
    public int getIlineDisplayInc() { return ilineDisplayInc; }
    public int getXlineDisplayInc() { return xlineDisplayInc; }
    public int getCorridorPercentage() { return corridorPercentage; }
    public boolean getShowAvgProfile() { return showAvgProfile; }
    public boolean getShowNeighborProfiles() { return showNeighborProfiles; }
    public boolean getShowLimitProfiles() { return showLimitProfiles; }
    public boolean getShowIntervalVelocities() { return showIntervalVelocities; }

    public void setShowLabels(boolean showLabels) { this.showLabels = showLabels; }
    public void setIlineDisplayInc(int displayInc) { ilineDisplayInc = displayInc; }
    public void setXlineDisplayInc(int displayInc) { xlineDisplayInc = displayInc; }
    public void setCorridorPercentage(int percentage) { corridorPercentage = percentage; }
    public void setShowAvgProfile(boolean show) { this.showAvgProfile = show; }
    public void setShowNeighborProfiles(boolean show) { this.showNeighborProfiles = show; }
    public void setShowLimitProfiles(boolean show) { this.showLimitProfiles = show; }
	public void setShowIntervalVelocities(boolean show) { this.showIntervalVelocities = show; }

}