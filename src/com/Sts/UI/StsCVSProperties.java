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

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

public class StsCVSProperties extends StsSeismicPanelProperties implements StsSerializable
{
    public int numberPanels = 11;
	public boolean ignoreSuperGather = false;
	public boolean autoStackUpdate = false;
    public boolean ignoreSuperGatherChanged = false;
    public int tracesPerPanel = 21;
    public int orientation = ORIENT_INLINE;
    public byte displayType = DISPLAY_WIGGLES;
    public StsWiggleDisplayProperties panelWiggleDisplayProperties;
   // public int panelWiggleOverlapPercent = 50;

    {
        //wiggleOverlapPercent = 30;
    }

    transient StsGroupBox groupBox = null;

	transient StsFloatFieldBean velocityStepBean;
    transient StsIntFieldBean numberPanelsBean;
    transient StsIntFieldBean tracesPerPanelBean;
	transient StsBooleanFieldBean ignoreSuperGatherBean;
	transient StsBooleanFieldBean autoStackUpdateBean;
    transient StsComboBoxFieldBean orientationBean;
    transient StsComboBoxFieldBean displayTypeBean;
    transient StsIntFieldBean wiggleOverlapPercentBean;

    transient StsSuperGatherProperties superProps;
	transient public boolean recompute = true;

	static private final String title = "CVS/VVS Properties";

    static public final String[] orientationStrings = new String[] { "Inline", "Crossline", "Angle up", "Angle down" };
    static public final int ORIENT_INLINE = 0;
    static public final int ORIENT_XLINE = 1;
    static public final int ORIENT_ANGLE_UP = 2;
    static public final int ORIENT_ANGLE_DOWN = 3;

    static public final String[] displayTypeStrings = new String[] { "Wiggles", "Var Area", "Var Density" };
    static public final int DISPLAY_WIGGLES = 0;
    static public final int DISPLAY_VAR_AREA = 1;
    static public final int DISPLAY_VAR_DENSITY = 2;

    public StsCVSProperties()
	{
	}

	public StsCVSProperties(String fieldName)
	{
		super(title, fieldName);
	}

	public StsCVSProperties(StsObject parentObject, StsCVSProperties defaultProperties, String fieldName)
	{
		super(parentObject, defaultProperties, title, fieldName);
	}

	public void initializeBeans()
	{
        ignoreSuperGatherBean = new StsBooleanFieldBean(this, "ignoreSuperGather", "Ignore Super Gather");
        autoStackUpdateBean = new StsBooleanFieldBean(this, "autoStackUpdate", "Auto-Update Stacks");
        tracesPerPanelBean = new StsIntFieldBean(this, "tracesPerPanel", 1, 999, "Number of Traces/Panel:", true);
        tracesPerPanelBean.setStep(2);
        StsPreStackLineSet lineSet = StsPreStackLineSetClass.currentProjectPreStackLineSet;
        if (lineSet == null) return;
        panelWiggleDisplayProperties = lineSet.getPanelWiggleDisplayProperties();
        wiggleOverlapPercentBean = new StsIntFieldBean(panelWiggleDisplayProperties, "wiggleOverlapPercent", 0, 10000, "Gain (Overlap Percent):", true); //sometimes higher gain needed SWC 6/8/09
        numberPanelsBean = new StsIntFieldBean(this, "numberPanels", 1, 100, "Number of Panels:", true);
        numberPanelsBean.setStep(2);
        displayTypeBean = new StsComboBoxFieldBean(this, "displayTypeString", "Display type:", displayTypeStrings);
        if(parentObject instanceof StsPreStackLineSet3d)
        {
            orientationBean = new StsComboBoxFieldBean(this, "orientationString", "Orientation", orientationStrings);
            propertyBeans = new StsFieldBean[] { orientationBean, tracesPerPanelBean, numberPanelsBean, wiggleOverlapPercentBean, displayTypeBean, ignoreSuperGatherBean, autoStackUpdateBean };
        }
        else
            propertyBeans = new StsFieldBean[] { tracesPerPanelBean, numberPanelsBean, wiggleOverlapPercentBean, displayTypeBean, ignoreSuperGatherBean, autoStackUpdateBean};
    }

    public int getNumberPanels() { return numberPanels; }
	public boolean getIgnoreSuperGather() { return ignoreSuperGather; }
    public int getTracesPerPanel() { return tracesPerPanel; }

    public void setNumberPanels(int nPanels)
    {
        if(nPanels <= 0) return;
        if(this.numberPanels == nPanels) return;
        numberPanels = nPanels;
        setRangeChanged(true);
    }

    public float getCvsVelocityStep(StsObject parentObject)
    {
        if(parentObject == null)
            parentObject = StsPreStackLineSetClass.currentProjectPreStackLineSet;
        velocityMin = ((StsPreStackLineSet) parentObject).semblanceRangeProperties.getVelocityMin() / 1000.0f;
        velocityMax = ((StsPreStackLineSet) parentObject).semblanceRangeProperties.getVelocityMax() / 1000.0f;
        return (velocityMax-velocityMin)/(float)(numberPanels-1);
    }

    public void setIgnoreSuperGather(boolean use)
    {
        ignoreSuperGather = use;
        //setRangeChanged(true);
        ignoreSuperGatherChanged = true;
        setChanged(true);
    }

    public void setTracesPerPanel(int ntraces)
    {
        tracesPerPanel = ntraces;
        setRangeChanged(true);
    }

    public void setOrientationString(String orientationString)
    {
        orientation = StsParameters.getStringMatchIndex(orientationStrings, orientationString);
        setChanged(true);
    }

    public String getOrientationString() { return orientationStrings[orientation]; }

    public int getOrientation() { return orientation; }


    public void setDisplayTypeString(String displayTypeString)
    {
        displayType = (byte)StsParameters.getStringMatchIndex(displayTypeStrings, displayTypeString);
        setWiggleParameters(displayType);
        setChanged(true);
    }

    private void setWiggleParameters(byte displayType2)
    {
        if(panelWiggleDisplayProperties == null) return;
        
        if (displayType2 == DISPLAY_WIGGLES)
        {
            panelWiggleDisplayProperties.setWiggleDrawLine(true);
            panelWiggleDisplayProperties.setWiggleDrawPoints(false);
        }
        else if (displayType2 == DISPLAY_VAR_AREA)
        {
            panelWiggleDisplayProperties.setWiggleDrawLine(false);
            panelWiggleDisplayProperties.setWiggleDrawPoints(false);
        }
        
    }

    public String getDisplayTypeString() { return displayTypeStrings[displayType]; }
    public byte getDisplayType() { return displayType; }

    public boolean ignoreSuperGatherChanged()
    {
        boolean changed = ignoreSuperGatherChanged;
        ignoreSuperGatherChanged = false;
        return changed;
    }
    
    public boolean getAutoStackUpdate()
    {
        return autoStackUpdate;
    }

    public void setAutoStackUpdate(boolean autoStackUpdate)
    {
        this.autoStackUpdate = autoStackUpdate;
        setChanged(true);
    }

}
