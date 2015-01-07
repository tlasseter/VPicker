package com.Sts.UI;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
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
public class StsSemblanceDisplayProperties extends StsPanelProperties implements StsSerializable
{
	public int mode = DISPLAY_MODE;
	public boolean showLabels = true;
	public boolean showAvgProfile = true;
	public boolean showNeighborProfiles = true;
	public boolean showLimitProfiles = true;
	public boolean showIntervalVelocities = true;
	public int ilineDisplayInc = 25;
	public int xlineDisplayInc = 25;
    public float corridorPercentage = 100;
//    public int traceThreshold = 10;
	public boolean displaySemblance = true;

	transient public StsModel model;
	transient protected StsPreStackLineSet3dClass seismicClass;

//	transient public StsJPanel propertyPanel = null; // When creating velocity profiles
	transient StsButtonFieldBean displayColorscaleButtonBean;
	transient StsButtonFieldBean editColorscaleButtonBean;

    transient StsButtonListFieldBean semblanceTypeBean;

    static public final int DISPLAY_MODE = StsPreStackLineSet3d.DISPLAY_MODE;
	static public final int EDIT_MODE = StsPreStackLineSet3d.EDIT_MODE;

	static public final String title = "Velocity Display Properties";
	static public final boolean immediateChange = true;

	transient public boolean redraw = false;

	public StsSemblanceDisplayProperties()
	{
	}

	public StsSemblanceDisplayProperties(int mode, String fieldName)
	{
		super(title, fieldName);
		this.mode = mode;
		ilineDisplayInc = 25;
		xlineDisplayInc = 25;
        corridorPercentage = StsPreStackLineSet.defaultCorridorPercentage;
//        traceThreshold = StsPreStackLineSet.defaultTraceThreshold;
	}

	public StsSemblanceDisplayProperties(StsObject parentObject, StsSemblanceDisplayProperties defaultProperties, int mode, String fieldName)
	{
        super(parentObject, title, fieldName);
        initializeDefaultProperties(defaultProperties);
        this.mode = mode;
		ilineDisplayInc = 25;
		xlineDisplayInc = 25;
        corridorPercentage = StsPreStackLineSet.defaultCorridorPercentage;
//        traceThreshold = StsPreStackLineSet.defaultTraceThreshold;
	}

    public void initializeDefaultProperties(Object defaultProperties)
    {
        if(defaultProperties != null)
            StsToolkit.copyObjectNonTransientFields(defaultProperties, this);
    }
	public void initializeBeans()
	{
		model = StsObject.getCurrentModel();
		seismicClass = (StsPreStackLineSet3dClass)model.getStsClass(StsPreStackLineSet3d.class);

		propertyBeans = new StsFieldBean[]
		{
			new StsStringFieldBean(this, "modeAsString", false, "Display or Edit Properties"),
			new StsIntFieldBean(this, "ilineDisplayInc", 0, 100, "Inline Display Distance (# of lines):", true),
			new StsIntFieldBean(this, "xlineDisplayInc", 0, 100, "Crossline Display Distance (# of gathers):", true),
            new StsBooleanFieldBean(this, "displaySemblance", "Show Semblance"),
			new StsBooleanFieldBean(this, "showLabels", "Show Point Labels"),
			new StsBooleanFieldBean(this, "showAvgProfile", "Show Avg Profile"),
			new StsBooleanFieldBean(this, "showNeighborProfiles", "Show Neighbor Profiles"),
			new StsBooleanFieldBean(this, "showLimitProfiles", "Show Limit Profiles"),
            new StsBooleanFieldBean(this, "showIntervalVelocities", "Show Interval Velocities"),
            new StsFloatFieldBean(this, "corridorPercentage", 10, 100, "Corridor Percentage:", true),
//            new StsIntFieldBean(this, "traceThreshold", 1, 25, "Trace Threshold:", true),
//			displayColorscaleButtonBean = new StsButtonFieldBean("Display colorscale", "Brings up colorscale for legend.", this, "displayColorscale"),
//			editColorscaleButtonBean = new StsButtonFieldBean("Edit colorscale", "Brings up colorscale for editing.", this, "editColorscale"),
		};
	}

	public int getMode() { return mode; }
    public String getModeAsString() { return StsPreStackLineSet3d.MODE_STRINGS[mode]; }
	public boolean getShowLabels() { return showLabels; }
	public int getIlineDisplayInc() { return ilineDisplayInc; }
	public int getXlineDisplayInc() { return xlineDisplayInc; }

    public float getCorridorPercentage() { return corridorPercentage; }
 //   public int getTraceThreshold() { return traceThreshold; }

	public boolean getDisplaySemblance() { return displaySemblance; }
	public boolean getShowAvgProfile() { return showAvgProfile; }
	public boolean getShowNeighborProfiles() { return showNeighborProfiles; }
	public boolean getShowLimitProfiles() { return showLimitProfiles; }
	public boolean getShowIntervalVelocities() { return showIntervalVelocities; }

    public void setMode(int mode) { this.mode = mode; }
	public void setShowLabels(boolean showLabels) { this.showLabels = showLabels; }
	public void setIlineDisplayInc(int displayInc) { ilineDisplayInc = displayInc; }
	public void setXlineDisplayInc(int displayInc) { xlineDisplayInc = displayInc; }

    public void setCorridorPercentage(float percentage) { corridorPercentage = percentage; }
 //   public void setTraceThreshold(int threshold) { traceThreshold = threshold; }

	public void setDisplaySemblance(boolean display) { displaySemblance = display; }
	public void setShowAvgProfile(boolean show) { this.showAvgProfile = show; }
	public void setShowNeighborProfiles(boolean show) { this.showNeighborProfiles = show; }
	public void setShowLimitProfiles(boolean show) { this.showLimitProfiles = show; }
	public void setShowIntervalVelocities(boolean show) { this.showIntervalVelocities = show; }
/*
    public void displayColorscale()
	{
		StsColorscalePanel colorPanel = new StsColorscalePanel(lineSetClass.semblanceColorList.colorscale, true);
		StsToolkit.createDialog(colorPanel, false);
	}

	public void editColorscale()
	{
		StsColorscalePanel colorPanel = new StsColorscalePanel(lineSetClass.semblanceColorList.colorscale, false);
		model.getActionManager().startAction(StsColorscaleAction.class, new Object[] { colorPanel, (Frame)model.win3d } );
//		redraw = true;
		model.win3dDisplayAll();
	}
*/
}
