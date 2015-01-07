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
public class StsSemblanceProperties extends StsPanelProperties implements StsSerializable
{
	public int mode = DISPLAY_MODE;
	public boolean showLabels = true;
	public boolean showAvgProfile = true;
	public boolean showNeighborProfiles = true;
	public boolean showLimitProfiles = true;
	public boolean showIntervalVelocities = true;
	public int ilineDisplayInc = 25;
	public int xlineDisplayInc = 25;
//    public int corridorPercentage = 100;
	public boolean displaySemblance = true;
	public byte order = ORDER_2ND;

	transient public StsModel model;
	transient protected StsPreStackLineSet3dClass seismicClass;

//	transient public StsJPanel propertyPanel = null; // When creating velocity profiles
//	transient StsButtonFieldBean displayColorscaleButtonBean;
//	transient StsButtonFieldBean editColorscaleButtonBean;

	static public final int DISPLAY_MODE = StsPreStackLineSet3d.DISPLAY_MODE;
	static public final int EDIT_MODE = StsPreStackLineSet3d.EDIT_MODE;

	/** order of NMO equation for semblance and gathers */
	static public final String ORDER_2ND_STRING = "Normal (2nd)";
	static public final String ORDER_4TH_STRING = "4th Order";
	static public final String ORDER_6TH_STRING = "6th Order";
	static public final String ORDER_6TH_OPT_STRING = "Opt 6th Order";
	static public final String[] orderStrings = new String[] { ORDER_2ND_STRING, ORDER_4TH_STRING, ORDER_6TH_STRING, ORDER_6TH_OPT_STRING };
	static public final byte ORDER_2ND = 0;
	static public final byte ORDER_4TH = 1;
	static public final byte ORDER_6TH = 2;
	static public final byte ORDER_6TH_OPT = 3;

	static public final String title = "Semblance Properties";
	static public final boolean immediateChange = true;

	transient public boolean redraw = false;

	public StsSemblanceProperties()
	{
	}

	public StsSemblanceProperties(int mode, String fieldName)
	{
		super(title, fieldName);
		this.mode = mode;
		ilineDisplayInc = 25;
		xlineDisplayInc = 25;
//        corridorPercentage = StsPreStackLineSet.defaultCorridorPercentage;
	}

	public StsSemblanceProperties(StsObject parentObject, StsSemblanceProperties defaultProperties, int mode, String fieldName)
	{
        super(parentObject, title, fieldName);
        initializeDefaultProperties(defaultProperties);
        this.mode = mode;
		ilineDisplayInc = 25;
		xlineDisplayInc = 25;
//        corridorPercentage = StsPreStackLineSet.defaultCorridorPercentage;
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
			new StsIntFieldBean(this, "ilineDisplayInc", 0, 100, "Inline Profile Display Distance:", true),
			new StsIntFieldBean(this, "xlineDisplayInc", 0, 100, "Crossline Profile Display Distance:", true),
			new StsBooleanFieldBean(this, "isDisplayData", "Show Semblance"),
			new StsBooleanFieldBean(this, "showLabels", "Show Point Labels"),
			new StsBooleanFieldBean(this, "showAvgProfile", "Show Avg Profile"),
			new StsBooleanFieldBean(this, "showNeighborProfiles", "Show Neighbor Profiles"),
			new StsBooleanFieldBean(this, "showLimitProfiles", "Show Limit Profiles"),
//            new StsIntFieldBean(this, "corridorPercentage", 10, 100, "Corridor Percentage:", true),
			new StsBooleanFieldBean(this, "showIntervalVelocities", "Show Interval Velocities"),
			new StsComboBoxFieldBean(this, "orderString", "Semblance Order", orderStrings),
//			displayColorscaleButtonBean = new StsButtonFieldBean("Display colorscale", "Brings up colorscale for legend.", this, "displayColorscale"),
//			editColorscaleButtonBean = new StsButtonFieldBean("Edit colorscale", "Brings up colorscale for editing.", this, "editColorscale"),
		};
	}

	public int getMode() { return mode; }
    public String getModeAsString() { return StsPreStackLineSet3d.MODE_STRINGS[mode]; }
	public boolean getShowLabels() { return showLabels; }
	public int getIlineDisplayInc() { return ilineDisplayInc; }
	public int getXlineDisplayInc() { return xlineDisplayInc; }
 //   public int getCorridorPercentage() { return corridorPercentage; }
	public boolean getDisplaySemblance() { return displaySemblance; }
	public boolean getShowAvgProfile() { return showAvgProfile; }
	public boolean getShowNeighborProfiles() { return showNeighborProfiles; }
	public boolean getShowLimitProfiles() { return showLimitProfiles; }
	public boolean getShowIntervalVelocities() { return showIntervalVelocities; }
	public String getOrderString() { return orderStrings[order]; }
	public byte getOrder() { return order; }

	public void setMode(int mode) { this.mode = mode; }
	public void setShowLabels(boolean showLabels) { this.showLabels = showLabels; }
	public void setIlineDisplayInc(int displayInc) { ilineDisplayInc = displayInc; }
	public void setXlineDisplayInc(int displayInc) { xlineDisplayInc = displayInc; }
 //   public void setCorridorPercentage(int percentage) { corridorPercentage = percentage; }
	public void setDisplaySemblance(boolean display) { displaySemblance = display; }
	public void setShowAvgProfile(boolean show) { this.showAvgProfile = show; }
	public void setShowNeighborProfiles(boolean show) { this.showNeighborProfiles = show; }
	public void setShowLimitProfiles(boolean show) { this.showLimitProfiles = show; }
	public void setShowIntervalVelocities(boolean show) { this.showIntervalVelocities = show; }
	public void setOrder(byte order) { this.order = order; }

	public void setOrderString(String orderString)
	{
		for(int n = 0; n < orderStrings.length; n++)
		{
			if(orderString == orderStrings[n])
			{
				this.order = (byte)n;
				return;
			}
		}
    }
/*
	public void displayColorscale()
	{
		StsColorscalePanel colorPanel = new StsColorscalePanel(lineSetClass.semblanceColorList.colorscale, true);
		StsToolkit.createDialog(colorPanel);
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
