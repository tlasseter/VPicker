package com.Sts.DBTypes;

import com.Sts.DB.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Nov 9, 2006
 * Time: 3:31:15 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsPreStackMicroseismicSetClass extends StsPreStackLineSetClass implements StsSerializable // , StsClassTextureDisplayable
{
	// public StsSemblanceDisplayProperties defaultSemblanceDisplayProperties;
	// public StsSemblanceDisplayProperties defaultSemblanceEditProperties;
	// public StsSemblanceComputeProperties defaultSemblanceComputeProperties;
	// public StsFilterProperties defaultFilterProperties;
	// public StsAGCPreStackProperties defaultAGCProperties;
    // public StsCVSProperties defaultCVSProperties;
	// public StsWiggleDisplayProperties defaultWiggleDisplayProperties = null;
	
    // private StsWiggleDisplayProperties wiggleDisplayProperties = null;
    
	// private String seismicSpectrumName = StsSpectrumClass.SPECTRUM_BWYR;
	// private String semblanceSpectrumName = StsSpectrumClass.SPECTRUM_SEMBLANCE;

	// protected boolean isPixelMode = false; // Blended Pixels or Nearest
	// protected boolean contourColors = true; // shader
	// protected boolean displayWiggles = true; // Display wiggles if data density allows - 2D Only
	// protected int wiggleToPixelRatio = 4;
	// private int displayPointSize = 10;
    // private boolean displayAxis = true;

	static public StsPreStackMicroseismicSet currentProjectPreStackMicroseismicSet = null;
    public StsPreStackMicroseismicSetClass()
	{
		super();
        userName = "Sets of Pre-Stack Microseismic Data";
	}

    public void initializePropertyFields()
	{
		propertyFields = new StsFieldBean[]
		{
			new StsButtonFieldBean("Wiggle Display Properties", "Edit wiggle display properties.", this, "displayWiggleProperties"),
			new StsButtonFieldBean("Semblance Display Properties", "Edit semblanceBytes display properties.", this, "displaySemblanceProperties"),
			new StsButtonFieldBean("AGC Properties", "Edit automatic gain control properties.", this, "displayAGCProperties")
		};
    }

    public void initializeDisplayFields()
	{
        displayFields = new StsFieldBean[]
		{
			new StsIntFieldBean(this, "displayPointSize", true, "Point Size:"),
			new StsBooleanFieldBean(this, "displayAxis", "Plot Axis"),
			new StsBooleanFieldBean(this, "isPixelMode", "Pixel Display Mode"),
			new StsBooleanFieldBean(this, "displayWiggles", "Wiggle Traces"),
			new StsBooleanFieldBean(this, "contourColors", "Contoured Prestack Colors"),
			new StsIntFieldBean(this, "wiggleToPixelRatio", 1, 100, "Wiggle to Pixel Ratio:")
		};
    }

    public void initializeDefaultFields()
    {
	    defaultFields =	new StsFieldBean[]
		{
			new StsComboBoxFieldBean(this, "seismicSpectrumName", "Seismic Spectrum:", StsSpectrumClass.cannedSpectrums),
			new StsComboBoxFieldBean(this, "semblanceSpectrumName", "Semblance Spectrum:", StsSpectrumClass.cannedSpectrums),
		};
	}
    
	public boolean projectInitialize(StsModel model)
	{
        defaultSemblanceDisplayProperties = new StsSemblanceDisplayProperties(StsPreStackLineSet3d.DISPLAY_MODE, "defaultSemblanceDisplayProperties");
		defaultSemblanceComputeProperties = new StsSemblanceComputeProperties("defaultSemblanceComputeProperties");
		defaultFilterProperties = new StsFilterProperties("defaultFilterProperties");
		defaultAGCProperties = new StsAGCPreStackProperties("defaultAGCProperties");
        defaultCVSProperties = new StsCVSProperties("defaultCVSProperties");
		defaultWiggleDisplayProperties = new StsWiggleDisplayProperties(this, "defaultWiggleDisplayProperties");
        return true;
	}
	
    public boolean dbInitialize()
    {
        /*
        if(wiggleDisplayProperties == null)
            wiggleDisplayProperties = new StsWiggleDisplayProperties(this, defaultWiggleDisplayProperties, "wiggleDisplayProperties");
        else
            wiggleDisplayProperties.setParentObject(this);
            */
        return true;
    }
    
    public void displayClass(StsGLPanel3d glPanel3d)
    {
        Iterator iter = getObjectIterator();
        while(iter.hasNext())
        {
        	StsPreStackMicroseismicSet microSet = (StsPreStackMicroseismicSet)iter.next();
            //microSet.display(glPanel3d);   // Not sure if we want to display something in 3d view for this type
        }
    }

	public boolean setCurrentObject(StsObject object)
	{
        StsPreStackMicroseismicSet newObject = (StsPreStackMicroseismicSet)object;
        if(newObject != null)
        {
            if(currentProjectPreStackMicroseismicSet == newObject) return false;
            currentProjectPreStackMicroseismicSet = newObject;
        }
        return super.setCurrentObject(object);
	}
    public StsObject getCurrentObject()
    {
        return currentProjectPreStackMicroseismicSet;
    }
	
	public String getSeismicSpectrumName()
	{
		return seismicSpectrumName;
	}

	public void setSeismicSpectrumName(String value)
	{
		if (this.seismicSpectrumName.equals(value))
			return;
		this.seismicSpectrumName = value;
	}

	public String getSemblanceSpectrumName()
	{
		return semblanceSpectrumName;
	}

	public void setSemblanceSpectrumName(String value)
	{
		if (this.semblanceSpectrumName.equals(value))
			return;
		this.semblanceSpectrumName = value;
	}
	public boolean getIsPixelMode()
	{
		return isPixelMode;
	}

	public void setIsPixelMode(boolean mode)
	{
		if (this.isPixelMode == mode)
			return;
		this.isPixelMode = mode;
		setDisplayField("isPixelMode", isPixelMode);
		currentModel.win3dDisplayAll();
	}

	public void setContourColors(boolean contourColors)
	{
		if (this.contourColors == contourColors) return;
		this.contourColors = contourColors;
		setDisplayField("contourColors", contourColors);
		list.forEach("semblanceColorListChanged");
		list.forEach("seismicColorListChanged");
	}
	
	public StsPreStackMicroseismicSet findMicroseismicSet(long eventTime)
	{
        Iterator iter = getObjectIterator();
        while(iter.hasNext())
        {
        	StsPreStackMicroseismicSet microSet = (StsPreStackMicroseismicSet)iter.next();
            if(microSet.isInTimeRange(eventTime))
            	return microSet;
        }		
        return null;
	}
	
	public boolean getContourColors()
	{
		return contourColors;
	}

	public boolean getDisplayWiggles()
	{
		return displayWiggles;
	}

	public void setDisplayWiggles(boolean mode)
	{
		if (this.displayWiggles == mode)
			return;
		this.displayWiggles = mode;
		setDisplayField("displayWiggles", displayWiggles);
		currentModel.win3dDisplayAll();
	}

	public int getWiggleToPixelRatio()
	{
		return wiggleToPixelRatio;
	}

	public void setWiggleToPixelRatio(int wtop)
	{
		if (this.wiggleToPixelRatio == wtop)
			return;
		this.wiggleToPixelRatio = wtop;
		setDisplayField("wiggleToPixelRatio", wiggleToPixelRatio);
		currentModel.win3dDisplayAll();
	}
	
    public StsWiggleDisplayProperties getWiggleDisplayProperties()
    {
        if(wiggleDisplayProperties == null)
            wiggleDisplayProperties = new StsWiggleDisplayProperties(this, defaultWiggleDisplayProperties, "wiggleDisplayProperties");
        return wiggleDisplayProperties;
    }
	public void displayWiggleProperties()
	{
		new StsOkApplyCancelDialog(currentModel.win3d, defaultWiggleDisplayProperties, "defaultWiggleDisplayProperties", false);
	}	
	public void displaySemblanceProperties()
	{
		new StsOkApplyCancelDialog(currentModel.win3d, this.defaultSemblanceDisplayProperties, "defaultDisplaySemblanceProperties", false);
	}
	public void displayAGCProperties()
	{
		new StsOkApplyCancelDialog(currentModel.win3d, defaultAGCProperties, "Edit AGC Properties", false);
	}
    public void displayFilterProperties()
    {
        new StsOkApplyCancelDialog(currentModel.win3d, defaultFilterProperties, "Edit Filter Properties", false);
	}
    public void displayCVSProperties()
    {
        new StsOkApplyCancelDialog(currentModel.win3d, defaultCVSProperties, "Edit CVS/VVS Properties", false);
	}

    public int getDisplayPointSize()
    {
        return displayPointSize;
    }

    public void setDisplayPointSize(int displayPointSize)
    {
        this.displayPointSize = displayPointSize;
    }

    public boolean isDisplayAxis()
    {
        return displayAxis;
    }

    public void setDisplayAxis(boolean displayAxis)
    {
        this.displayAxis = displayAxis;
    }
}
