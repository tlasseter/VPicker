package com.Sts.DBTypes;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.DB.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

import javax.swing.*;
import java.awt.*;

public class StsVspClass extends StsClass implements StsSerializable, StsClassDisplayable //, StsEfficientRepaintable
{
	protected boolean displayAxis = true;
	protected boolean displayVAR = true;
	protected boolean isPixelMode = false; // Blended Pixels or Nearest
	//protected boolean displayWiggles = false; // Display wiggles if data density allows - 2D Only
	protected int wiggleToPixelRatio = 4;
	protected boolean displayVSPs = true;
	protected boolean contourColors = true;
	protected float tracesPerInch = 25.f;
	public StsWiggleDisplayProperties defaultWiggleProperties = null;
	public StsFilterProperties defaultFilterProperties = null;
	public StsAGCProperties defaultAGCProperties = null;
    public StsDatumProperties defaultDatumProperties = null;
	private String seismicSpectrumName = StsSpectrumClass.SPECTRUM_RWB;

   /** display properties for all traces associated with this volume */
    private StsWiggleDisplayProperties wiggleDisplayProperties = null;

    static public final String VSP_RAW_STRING = "Raw";
	static public final String VSP_FLDSTACK_STRING = "FieldStack";
	static public final String VSP_PROCESSED_STRING = "Processed";
	static public final String VSP_MIGRATED_STRING = "Migrated";
	static public final String VSP_SYN_STRING = "Synthetic";
	public static String[] vspTypeStrings = { VSP_RAW_STRING, VSP_FLDSTACK_STRING, VSP_PROCESSED_STRING, VSP_MIGRATED_STRING, VSP_SYN_STRING };

	static public final byte VSP_RAW = 0;
	static public final byte VSP_FLDSTACK = 1;
	static public final byte VSP_PROCESSED = 2;
	static public final byte VSP_MIGRATED = 3;
	static public final byte VSP_SYNTHETIC = 4;

    protected boolean enableTime = true;

    transient StsSubType[] subTypes = null;

    // TODO should move properties constructor to projectInitialize; default constructor should be a null constructor for DB
    public StsVspClass()
	{
        userName = "Vertical Seismic Profiles";
    }

    public boolean projectInitialize(StsModel model)
	{
		defaultWiggleProperties = new StsWiggleDisplayProperties(this, "defaultWiggleProperties");
		defaultFilterProperties = new StsFilterProperties("defaultFilterProperties");
		defaultAGCProperties = new StsAGCProperties("defaultAGCProperties");
        defaultDatumProperties = new StsDatumProperties("defaultDatumProperties");
        wiggleDisplayProperties = new StsWiggleDisplayProperties(this, defaultWiggleProperties, "wiggleDisplayProperties");
        initializeSubClasses();
        return true;
    }
    
    private void initializeSubClasses()
    {
	    subTypes = new StsSubType[vspTypeStrings.length];
		for(int n = 0; n < vspTypeStrings.length; n++)
			subTypes[n] = new StsSubType(this, vspTypeStrings[n], (byte)n);
    }

    public boolean dbInitialize()
    {
        if(wiggleDisplayProperties == null)
            wiggleDisplayProperties = new StsWiggleDisplayProperties(this, defaultWiggleProperties, "wiggleDisplayProperties");
        else
            wiggleDisplayProperties.setParentObject(this);
        return true;
    }

    public void initializeDisplayFields()
	{

		displayFields = new StsFieldBean[]
			{
			new StsBooleanFieldBean(this, "displayVSPs", "Display on Well"),
            new StsBooleanFieldBean(this, "enableTime", "Enable Time"),
			new StsBooleanFieldBean(this, "displayAxis", "Plot Axis"),
			new StsBooleanFieldBean(this, "isPixelMode", "Pixel Display Mode"),
			new StsBooleanFieldBean(this, "displayVAR", "Variable Area"),
			new StsBooleanFieldBean(this, "contourColors", "Contoured Seismic Colors"),
			new StsIntFieldBean(this, "wiggleToPixelRatio", 1, 100, "Minimum Pixels per Wiggle"),
			//new StsFloatFieldBean(this, "tracesPerInch",1,100,"Traces per Inch:")
		};

	}

	public void initializeDefaultFields()
	{
        defaultFields = new StsFieldBean[]
			{
			new StsComboBoxFieldBean(this, "seismicSpectrumName", "Seismic Spectrum:", StsSpectrumClass.cannedSpectrums),
			new StsButtonFieldBean("Wiggle Display Properties", "Default wiggle display properties.", this, "displayWiggleProperties"),

		};
	}

	public StsOkApplyCancelDialog displayWigglePropertiesDialog(String title, Frame frame)
	{

		return new StsOkApplyCancelDialog(frame, defaultWiggleProperties, title, false);
	}

	public void displayWiggleProperties()
    {
	   defaultWiggleProperties.displayWiggleProperties("Edit Default VSP wiggle properties");
    }

    public void displayAGCProperties()
    {
        new StsOkApplyCancelDialog(currentModel.win3d, defaultAGCProperties, "Edit AGC Display Properties", false);
    }

    public void displayDatumProperties()
    {
            new StsOkApplyCancelDialog(currentModel.win3d, defaultDatumProperties, "Edit Datum Display Properties", false);
	}
	public void selected(StsVsp vsp)
	{
		super.selected(vsp);
		setCurrentObject(vsp);
	}

	/**
	 * Set the pixel mode for VSP displays. When set to off, pixels will be interpolated between samples.
	 * When set to on, there will be hard edges between samples.
	 * @param isPixelMode boolean
	 */
	public void setIsPixelMode(boolean isPixelMode)
	{
		if(this.isPixelMode == isPixelMode)return;
		this.isPixelMode = isPixelMode;
		//setDisplayField("isPixelMode", isPixelMode);
		forEach("setIsPixelMode",new Boolean(isPixelMode));
		//currentModel.win3dDisplayAll();
	}

	public boolean getIsPixelMode()
	{
		return isPixelMode;
	}

    public void setEnableTime(boolean enable)
    {
        if(this.enableTime == enable) return;
        this.enableTime = enable;
//        setDisplayField("enableTime", enableTime);
        currentModel.win3dDisplayAll();
    }
    public boolean getEnableTime() {	return enableTime; }
	/**
	 * Set the display axis toggle
	 */
	public void setDisplayAxis(boolean displayAxis)
	{
		if(this.displayAxis == displayAxis)return;
		this.displayAxis= displayAxis;
		//setDisplayField("displayAxis", displayAxis);
		forEach("setDisplayAxis",new Boolean(displayAxis));
		//currentModel.win3dDisplayAll();
	}

	public boolean getDisplayAxis()
	{
		return displayAxis;
	}

	public boolean getContourColors()
	{
		return contourColors;
	}
	/** display shader */
	public void setContourColors(boolean val)
	{
		if(this.contourColors == val)return;
		this.contourColors = val;
		forEach("setContourColors",new Boolean(contourColors));
	}

	/** display variable area toggle */
	public void setDisplayVAR(boolean displayVAR)
	{
		if(this.displayVAR == displayVAR)return;
		this.displayVAR = displayVAR;
		//setDisplayField("displayVAR", displayVAR);
		forEach("setDisplayVAR",new Boolean(displayVAR));
		//currentModel.win3dDisplayAll();
	}

	public boolean getDisplayVAR()
	{
		return displayVAR;
	}

	/**
	 * If the display wiggle flag is set to on, wiggle traces will be plotted on top of the
	 * texture maps and when they are isVisible is based on the wiggle to pixel ratio. The wiggle
	 * to pixel ratio is the number of pixels between traces before wiggles will be plotted. This is
	 * no avoid ridiculous plots where the number of traces plotted exceeds number of pixles resulting in
	 * solid black plots.
	 * @param wiggleToPixelRatio int
	 */
	public void setWiggleToPixelRatio(int wiggleToPixelRatio)
	{
		if(this.wiggleToPixelRatio == wiggleToPixelRatio)return;
		this.wiggleToPixelRatio = wiggleToPixelRatio;
		forEach("setWiggleToPixelRatio",new Integer(wiggleToPixelRatio));
		currentModel.win3dDisplayAll();
	}

	public int getWiggleToPixelRatio()
	{
		return wiggleToPixelRatio;
	}

/*
	public void setTracesPerInch(float tracesPerInch)
	{
		if(this.tracesPerInch == tracesPerInch)return;
		this.tracesPerInch = tracesPerInch;
		forEach("setTracesPerInch",new Float(tracesPerInch));
		//setDisplayField("tracesPerInch", tracesPerInch);
		//currentModel.win3dDisplayAll();
	}

	public float getTracesPerInch()
	{
		return tracesPerInch;
	}
*/
   public void setTracesPerInch(float tracesPerInch)
   {
	 defaultWiggleProperties.setTracesPerInch(tracesPerInch);
   }

   public float getTracesPerInch()
   {
	   return defaultWiggleProperties.getTracesPerInch();
   }
   public void setInchesPerSecond(float inchesPerSecond)
   {
     defaultWiggleProperties.setInchesPerSecond(inchesPerSecond);
   }

   public float getInchesPerSecond()
   {
       return defaultWiggleProperties.getInchesPerSecond();
   }

	public void setSeismicSpectrumName(String seismicSpectrumName)
	{
		if(this.seismicSpectrumName == seismicSpectrumName)return;
		this.seismicSpectrumName = seismicSpectrumName;
		//setDisplayField("seismicSpectrumName", seismicSpectrumName);
		//currentModel.win3dDisplayAll();
	}

	public String getSeismicSpectrumName()
	{
		return seismicSpectrumName;
	}

    // efficientRepaintable form wiggle props -- wiggle props have changed!
	// propagate down
/*
	public void repaint()
	{
		for(int n = 0; n < getSize(); n++)
		{
			StsVsp vsp = (StsVsp)getElement(n);
			vsp.repaint();
		}

    }
*/
	public void displayTimeClass(StsGLPanel3d glPanel3d, long time)
	{
		for(int n = 0; n < getSize(); n++)
		{
			StsVsp vsp = (StsVsp)getElement(n);
            if((enableTime && vsp.isAlive(time)) || (!enableTime))
                vsp.display(glPanel3d);
		}
	}

	public boolean getDisplayVSPs()
	{
		return displayVSPs;
	}

	public void setDisplayVSPs(boolean displayVSPs)
	{
		if(this.displayVSPs == displayVSPs)return;
		this.displayVSPs = displayVSPs;
//		setDisplayField("displayVSPs", displayVSPs);
		forEach("setDisplayVSPs",new Boolean(displayVSPs));
		currentModel.win3dDisplayAll();
	}


	public StsSubType[] getSubTypes()
	{
		return subTypes;
	}

	static public byte getTypeForName(String typename)
	{
		for(byte n = 0; n < vspTypeStrings.length; n++)
			if(vspTypeStrings[n] == typename) return n;
		return (byte)0;
	}

    public boolean hasVSP(StsWell well)
    {
        for(int n = 0; n < getSize(); n++)
        {
            StsVsp vsp = (StsVsp)getElement(n);
            if(vsp.getWell() == well)
                return true;
		}
        return false;
    }

    public void createPopupMenu(StsWell well, JMenu menu)
    {
        JMenu subMenu = null;
        StsMenuItem menuItem;

        int nVsps = getSize();
        for(int n = 0; n < nVsps; n++)
        {
            StsVsp vsp = (StsVsp)getElement(n);
            if(vsp.getWell() == well)
            {
                subMenu = new JMenu(vsp.getName());
                menuItem = new StsMenuItem();
                menuItem.setMenuActionListener("Export...", vsp, "export", null);
                subMenu.add(menuItem);
            }
            if(subMenu != null)
                menu.add(subMenu);
        }
    }

    public StsWiggleDisplayProperties getWiggleDisplayProperties()
    {
        return wiggleDisplayProperties;
    }
}