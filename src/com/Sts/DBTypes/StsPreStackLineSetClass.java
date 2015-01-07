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
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;
import com.Sts.Utilities.*;

public class StsPreStackLineSetClass extends StsClass implements StsSerializable
{
	protected boolean displayGridLines = false;
	protected boolean displayAxis = true;
	protected boolean applySurfacesSubVolume = false;
	public String seismicSpectrumName = StsSpectrumClass.SPECTRUM_RWB;
	public String semblanceSpectrumName = StsSpectrumClass.SPECTRUM_SEMBLANCE;
    private String stacksSpectrumName = StsSpectrumClass.SPECTRUM_STACKS;
    private boolean showVelStat = true;
	protected int displayPointSize = 6;
	//protected boolean displayVelocity = true;
	private byte displayType = DISPLAY_VELOCITY;

	protected float stackNeighborRadius = 1000.0f;
	protected int gatherIntervalToStack = 1;
	protected float sampleIntervalToStack = 2;
	protected float maximumZToStack = 5000;
	protected float maximumOffsetToStack = 1000;

	/** skip gathers that have already been picked or unpicked or go through all */
	protected byte profileOption = EVERY_PROFILE;

	/** type of stacking to do */
	protected byte stackOption = STACK_LINES;

	protected boolean isPixelMode = false; // Blended Pixels or Nearest
	protected boolean contourColors = true; // shader
	protected boolean displayWiggles = true; // Display wiggles if data density allows - 2D Only
	protected int wiggleToPixelRatio = 4;

	public StsSemblanceDisplayProperties defaultSemblanceDisplayProperties;
	public StsSemblanceDisplayProperties defaultSemblanceEditProperties;
	public StsSemblanceComputeProperties defaultSemblanceComputeProperties;
//	public StsSemblanceRangeProperties defaultSemblanceRangeProperties;
	public StsFilterProperties defaultFilterProperties;
	public StsAGCPreStackProperties defaultAGCProperties;
    public StsDatumProperties defaultDatumProperties;
    public StsSuperGatherProperties defaultSuperGatherProperties;
    public StsCVSProperties defaultCVSProperties;
	public StsWiggleDisplayProperties defaultWiggleDisplayProperties = null;

    /** display properties for all traces associated with this volume */
    public StsWiggleDisplayProperties wiggleDisplayProperties = null;
    public StsWiggleDisplayProperties panelWiggleDisplayProperties;

//	transient private StsFloatFieldBean tpiBean = null;
//	transient public StsColorList semblanceColorList = null;

	private boolean flatten = false;

	static public StsPreStackLineSet currentProjectPreStackLineSet = null;

	static public final String EVERY_PROFILE_STRING = "All";
	static public final String PICKED_PROFILES_STRING = "Picked";
	static public final String UNPICKED_PROFILES_STRING = "Unpicked";
	static public final byte EVERY_PROFILE = 0;
	static public final byte PICKED_PROFILES = 1;
	static public final byte UNPICKED_PROFILES = 2;
	static public final String[] PROFILE_STRINGS = new String[] { EVERY_PROFILE_STRING, PICKED_PROFILES_STRING, UNPICKED_PROFILES_STRING };

    static public final String SINGLE_GATHER_STRING = "Single Gather";
    static public final String SUPER_CROSS_STRING = "Super Gather - Cross";
	static public final String SUPER_RECT_STRING = "Super Gather - Rectangular";
    static public final String SUPER_INLINE_STRING = "Super Gather - Inline";
    static public final String SUPER_XLINE_STRING = "Super Gather - Xline";
    static public final byte SUPER_SINGLE = 0;
    static public final byte SUPER_CROSS = 1;
	static public final byte SUPER_RECT = 2;
    static public final byte SUPER_INLINE = 3;
    static public final byte SUPER_XLINE = 4;
    static public final String[] GATHER_TYPE_STRINGS = new String[] { SINGLE_GATHER_STRING, SUPER_CROSS_STRING, SUPER_RECT_STRING, SUPER_INLINE_STRING, SUPER_XLINE_STRING };

	static public final String STACK_NONE_STRING = "None";
	static public final String STACK_LINES_STRING = "Lines";
	static public final String STACK_NEIGHBORS_STRING = "Neighbors";
//	static public final String STACK_VOLUME_STRING = "PostStack3d";
	static public final byte STACK_NONE = 0;
	static public final byte STACK_NEIGHBORS = 1;
	static public final byte STACK_LINES = 2;
	static public final byte STACK_VOLUME = 3;
	static public final String[] STACK_STRINGS = new String[] { STACK_NONE_STRING, STACK_NEIGHBORS_STRING, STACK_LINES_STRING, /*STACK_VOLUME_STRING*/ };

    static public final byte DISPLAY_VELOCITY = 0;
	static public final byte DISPLAY_STACKED = 1;
	static public final byte DISPLAY_SEMBLANCE = 2;
    static public final byte DISPLAY_NONE = 3;
    static public final byte DISPLAY_ATTRIBUTE = -1;

    static public final int nVolumeTypes = 3;

	static public final String VELOCITY_STRING = "Velocity";
	static public final String STACKED_STRING = "Stacked";
	static public final String SEMBLANCE_STRING = "Semblance";
    static public final String NONE_STRING = "None";
    static public final String[] DISPLAY_TYPE_STRINGS = new String[] { VELOCITY_STRING, STACKED_STRING, SEMBLANCE_STRING, NONE_STRING };

//	private StsColorscale semblanceColorscale = null;
//	private StsGLColorscaleDisplayList semblanceColorscaleDisplayList = null;

	public StsPreStackLineSetClass()
	{
	}

	public boolean projectInitialize(StsModel model)
	{
        defaultSemblanceDisplayProperties = new StsSemblanceDisplayProperties(StsPreStackLineSet3d.DISPLAY_MODE, "defaultSemblanceDisplayProperties");
		defaultSemblanceEditProperties = new StsSemblanceDisplayProperties(StsPreStackLineSet3d.EDIT_MODE, "defaultSemblanceDisplayProperties");
		defaultSemblanceComputeProperties = new StsSemblanceComputeProperties("defaultSemblanceComputeProperties");

//		defaultSemblanceRangeProperties = new StsSemblanceRangeProperties(model, this, "defaultSemblanceRangeProperties");
		defaultFilterProperties = new StsFilterProperties("defaultFilterProperties");
		defaultAGCProperties = new StsAGCPreStackProperties("defaultAGCProperties");
        defaultCVSProperties = new StsCVSProperties("defaultCVSProperties");
        defaultDatumProperties = new StsDatumProperties("defaultDatumProperties");
        defaultSuperGatherProperties = new StsSuperGatherProperties(this, "defaultSuperGatherProperties");
		defaultWiggleDisplayProperties = new StsWiggleDisplayProperties(this, "defaultWiggleDisplayProperties");
        return true;
	}

    public boolean dbInitialize()
    {
        if(wiggleDisplayProperties == null)
            wiggleDisplayProperties = new StsWiggleDisplayProperties(this, defaultWiggleDisplayProperties, "wiggleDisplayProperties");
        else
            wiggleDisplayProperties.setParentObject(this);
        return true;
    }

    public void finalInitialize()
	{
        if(currentProjectPreStackLineSet != null)
        {
            StsVelocityAnalysisToolbar.checkAddToolbar(currentModel, currentProjectPreStackLineSet, true);
            currentProjectPreStackLineSet.initializeVelStatusPoints();
        }
	}

	public void selected(StsPreStackLineSet lineSet)
	{
		super.selected(lineSet);
		setCurrentObject(lineSet);
	}

	public boolean setCurrentObject(StsObject object)
	{
        StsPreStackLineSet newObject = (StsPreStackLineSet)object;
        if(newObject != null)
        {
            if(currentProjectPreStackLineSet == newObject) return false;
            if (currentProjectPreStackLineSet != null) currentProjectPreStackLineSet.setShowVelStatusPoints(false);
            currentProjectPreStackLineSet = newObject;
            currentProjectPreStackLineSet.setShowVelStatusPoints(true);
            StsVelocityAnalysisToolbar.checkAddToolbar(currentModel, currentProjectPreStackLineSet, true);
        }
        return super.setCurrentObject(object);
	}

    public StsObject getCurrentObject()
    {
        return currentProjectPreStackLineSet;
    }
/*
    public void setViewObject(StsObject object)
    {
        currentModel.viewObjectChangedAndRepaint(object);
    }
*/
    /*
	 public boolean getDisplayGridLines() { return displayGridLines; }
	 public void setDisplayGridLines(boolean isLines)
	 {
	  if(this.displayGridLines == isLines) return;
	  this.displayGridLines = isLines;
	  setDisplayField("displayGridLines", displayGridLines);
	  currentModel.win3dDisplayAll();
	 }
	 */
	public boolean getDisplayAxis()
	{
		return displayAxis;
	}

	public void setDisplayAxis(boolean b)
	{
		if (this.displayAxis == b)
			return;
		this.displayAxis = b;
//		setDisplayField("displayAxis", displayAxis);
		currentModel.win3dDisplayAll();
	}

	public boolean getApplySurfacesSubVolume()
	{
		return applySurfacesSubVolume;
	}

	public void setApplySurfacesSubVolume(boolean apply)
	{
		if (this.applySurfacesSubVolume == apply)
			return;
		this.applySurfacesSubVolume = apply;
//		setDisplayField("applySurfacesSubVolume", applySurfacesSubVolume);
		currentModel.win3dDisplayAll();
	}

	public void setDisplayPointSize(int size)
	{
		if (this.displayPointSize == size)
			return;
		this.displayPointSize = size;
//		setDisplayField("displayPointSize", displayPointSize);
		currentModel.win3dDisplayAll();
	}

	public int getDisplayPointSize()
	{
		return displayPointSize;
	}

	// if currentSeismicVolume has been deleted, reset currentSeismicVolume to first available

	public StsPreStackLineSet getCurrentProjectLineSet()
	{
		return currentProjectPreStackLineSet;
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
    public void displayDatumProperties()
    {
            new StsOkApplyCancelDialog(currentModel.win3d, defaultDatumProperties, "Edit Datum Properties", false);
	}
    public void displaySuperGatherProperties()
    {
            new StsOkApplyCancelDialog(currentModel.win3d, defaultSuperGatherProperties, "Edit Super Gather Properties", false);
	}
	/*
	 public boolean setCurrentSeismicVolumeName(String name)
	 {
	  StsPreStackLineSet newSeismicVolume = (StsPreStackLineSet)getObjectWithName(name);
	  return setCurrentObject(newSeismicVolume);
	 }

	 public StsPreStackLineSet getSeismicVolumeDisplayableOnSurface()
	 {
	  if(currentObject == null) return null;
	  if(!isVisibleOnSurface || !isVisible) return null;
	  return (StsPreStackLineSet)currentObject;
	 }
	 */

	/*
	 public void displayOnCursor(StsCursor3d.StsCursor3dTexture cursorSection, StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, boolean is3d)
               {
                if(currentObject == null || !((StsPreStackLineSet)currentObject).getIsVisibleOnCursor()) return;
                StsPreStackSeismicCursorSection seismicCursorSection = cursorSection.hasCursor3dDisplayable(StsPreStackSeismicCursorSection.class);
                if(seismicCursorSection == null)
                {
                 seismicCursorSection = new StsPreStackSeismicCursorSection((StsPreStackLineSet)currentObject, glPanel3d, dirNo, dirCoordinate);
                 cursorSection.addDisplayable(seismicCursorSection);
                }
                seismicCursorSection.display(this, glPanel3d, is3d);
               }
               */
	/*
	 public void displayOnCursor2d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate)
	 {
	  if(currentObject == null || !((StsPreStackLineSet)currentObject).getIsVisibleOnCursor()) return;
	  GL gl = glPanel3d.getGL();
	  GLU glu = glPanel3d.getGLU();
//        currentSeismicVolume.displayOnCursor2d(dirNo, dirCoordinate, gl, glu);
	 }
	 */

	/*
	 public boolean getApplyCrop() { return applyCrop; }
	 public void setApplyCrop(boolean crop)
	 {
	  if(crop == applyCrop) return;
	  applyCrop = crop;
	  currentModel.getProject().getCropVolume().valueChanged();
	  currentModel.clearDataDisplays();
	  currentModel.win3dDisplayAll();
	 }
	 */
	public boolean getIsPixelMode()
	{
		return isPixelMode;
	}

	public void setIsPixelMode(boolean mode)
	{
		if (this.isPixelMode == mode)
			return;
		this.isPixelMode = mode;
//		setDisplayField("isPixelMode", isPixelMode);
		currentModel.win3dDisplayAll();
	}

	public void setContourColors(boolean contourColors)
	{
		if (this.contourColors == contourColors) return;
		this.contourColors = contourColors;
//		setDisplayField("contourColors", contourColors);
		list.forEach("semblanceColorListChanged");
		list.forEach("seismicColorListChanged");
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
//		setDisplayField("displayWiggles", displayWiggles);
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
//		setDisplayField("wiggleToPixelRatio", wiggleToPixelRatio);
		currentModel.win3dDisplayAll();
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

    public String getStacksSpectrumName()
    {
        return stacksSpectrumName;
    }

    public void setStacksSpectrumName(String value)
    {
        if (this.stacksSpectrumName.equals(value))
            return;
        this.stacksSpectrumName = value;
    }

	public boolean getShowVelStat()
	{
		return showVelStat;
	}

	public void setShowVelStat(boolean val)
	{
		showVelStat = val;
        if(currentProjectPreStackLineSet == null) return;
//        setDisplayField("showVelStat", showVelStat);

        currentProjectPreStackLineSet.setShowVelStatusPoints(val);
		currentModel.repaintViews(StsView3d.class);
    }

	public int getGatherIntervalToStack()
	{
		return gatherIntervalToStack;
	}

	public void setGatherIntervalToStack(int val)
	{
		if (this.gatherIntervalToStack == val)
			return;
		this.gatherIntervalToStack = val;
	}

	public float getSampleIntervalToStack()
	{
		return sampleIntervalToStack;
	}

	public void setSampleIntervalToStack(float val)
	{
		if (this.sampleIntervalToStack == val)
			return;
		this.sampleIntervalToStack = val;
	}

	public float getMaxZToStack()
	{
		return maximumZToStack;
	}

	public void setMaxZToStack(float val)
	{
		if (this.maximumZToStack == val)
			return;
		this.maximumZToStack = val;
	}

	public float getMaxOffsetToStack()
	{
		return maximumOffsetToStack;
	}

	public void setMaxOffsetToStack(float val)
	{
		if (this.maximumOffsetToStack == val)
			return;
		this.maximumOffsetToStack = val;
//		setDisplayField("maximumOffsetToStack", maximumOffsetToStack);
	}

	/*public float getTracesPerInch()
	  {
	 return tracesPerInch;
	  }

	  public void setTracesPerInch(float val)
	  {
	 if(this.tracesPerInch == val)
			   return;
	 this.tracesPerInch = val;
	 setDisplayField("tracesPerInch", tracesPerInch);
	 currentModel.win3dDisplayAll();
	  }

	  public void setTracesPerInch(float val, boolean update)
	 {
	  if(this.tracesPerInch == val)
	   return;
	  this.tracesPerInch = val;
	  setDisplayField("tracesPerInch", tracesPerInch);
	  if (update && tpiBean != null)
	   tpiBean.setValue(tracesPerInch);
	  currentModel.win3dDisplayAll();
	  }
	 */

// none of methods below are currently being used
/*
    public void setTracesPerInch(float val)
	{
		defaultWiggleDisplayProperties.setTracesPerInch(val);
		currentModel.win3dDisplayAll();
	}

	public void setTracesPerInch(float val, boolean update)
	{
		defaultWiggleDisplayProperties.setTracesPerInch(val);
		currentModel.win3dDisplayAll();
	}

	public float getTracesPerInch()
	{
		return defaultWiggleDisplayProperties.getTracesPerInch();
    }

    public int getNumPanels()
    {
        return defaultCVSProperties.getNumberPanels();
    }

    public float getMinCVSVelocity()
    {
        return defaultCVSProperties.getVelocityMin();
    }

    public float getMaxCVSVelocity()
    {
        return defaultCVSProperties.getVelocityMax();
    }

    public int getTracesPerPanel()
    {
        return defaultCVSProperties.getTracesPerPanel();
    }

	public float getInchesPerSecond()
	{
		return defaultWiggleDisplayProperties.getInchesPerSecond();
	}

	public void setInchesPerSecond(float val)
	{
		defaultWiggleDisplayProperties.setInchesPerSecond(val);
		currentModel.win3dDisplayAll();
	}

	public void setInchesPerSecond(float val, boolean update)
	{
		defaultWiggleDisplayProperties.setInchesPerSecond(val);
		currentModel.win3dDisplayAll();
	}
*/
	public void close()
	{
		list.forEach("close");
	}

	public boolean drawLast()
	{
		return false;
	}

	public boolean setFlatten()
	{
        if(!canFlatten()) return false;
        if (flatten) return true;
        flatten = true;
        currentModel.repaintViews(StsViewGather.class);
        return true;
    }

    private boolean canFlatten()
    {
        if(currentProjectPreStackLineSet == null) return false;
        return currentProjectPreStackLineSet.canFlatten();
    }

    public boolean setUnFlatten()
	{
        if(!canUnflatten()) return false;
        if (!flatten) return true;
		flatten = false;
		currentModel.repaintViews(StsViewGather.class);
        return true;
    }

    private boolean canUnflatten()
    {
        if(currentProjectPreStackLineSet == null) return false;
        return currentProjectPreStackLineSet.canUnflatten();
    }

	public boolean getFlatten()
	{
		return flatten;
	}

	public byte getProfileOption()
	{
		return profileOption;
	}

	public String getProfileOptionString()
	{
		return PROFILE_STRINGS[profileOption];
	}

	public void setProfileOptionString(String option)
	{
		for (byte i = 0; i < PROFILE_STRINGS.length; i++)
		{
			if (PROFILE_STRINGS[i] == option)
			{
				profileOption = i;
				return;
			}
		}
	}

	public byte getStackOption()
	{
		return stackOption;
	}

	public void setStackOption(byte option)
	{
		this.stackOption = option;
		int nSeismicVolumes = this.getSize();
		for (int n = 0; n < nSeismicVolumes; n++)
		{
			StsPreStackLineSet lineSet = (StsPreStackLineSet) getElement(n);
			lineSet.setStackOption(option);
		}
		currentModel.viewObjectRepaint(this, currentProjectPreStackLineSet);
	}

	public String getStackOptionString()
	{
		return STACK_STRINGS[stackOption];
	}

	public void setStackOptionString(String option)
	{
		for (byte i = 0; i < STACK_STRINGS.length; i++)
		{
			if (STACK_STRINGS[i].equals(option))
			{
				setStackOption(i);
				return;
			}
		}
	}

	public float getStackNeighborRadius()
	{
		return stackNeighborRadius;
	}

	public void setStackNeighborRadius(float radius)
	{
		stackNeighborRadius = radius;
	}

	/*
	 public void displayVelocity()
	 {
	  if(displayVelocity) return;
	  displayVelocity = true;
	  currentModel.win3dDisplayAll();
	 }
	 */

	public byte getDisplayType()
	{
		return displayType;
	}

	public boolean setDisplayType(final byte displayType)
	{
		if (this.displayType == displayType) return false;
		if (currentModel.win3d == null) return false;
		final StsComboBoxFieldBean displayTypeCombo = (StsComboBoxFieldBean)currentModel.win3d.getToolbarComponentNamed(StsVelocityAnalysisToolbar.NAME, StsVelocityAnalysisToolbar.DISPLAY_TYPE);
		if (displayTypeCombo != null)
		{
            Runnable runnable = new Runnable()
            {
 				public void run()
                {
                    displayTypeCombo.setSelectedIndex(displayType);
                }
            };
            StsToolkit.runLaterOnEventThread(runnable);
		}
		this.displayType = displayType;
		currentModel.viewObjectRepaint(this, currentProjectPreStackLineSet);
		return true;
	}

	public void computeStackedVolume()
	{
        if(currentProjectPreStackLineSet == null) return;
        if(!(currentProjectPreStackLineSet instanceof StsPreStackLineSet3d)) return;
		currentProjectPreStackLineSet.computePreStackVolume(DISPLAY_STACKED);
	}

	public void computeSemblanceVolume()
	{
        if(currentProjectPreStackLineSet == null) return;
        if(!(currentProjectPreStackLineSet instanceof StsPreStackLineSet3d)) return;
		currentProjectPreStackLineSet.computePreStackVolume(DISPLAY_SEMBLANCE);
	}

	public String getDisplayTypeString()
	{
		return DISPLAY_TYPE_STRINGS[displayType];
	}

	public void setDisplayTypeString(String displayType)
	{
		for (byte i = 0; i < DISPLAY_TYPE_STRINGS.length; i++)
		{
			if (DISPLAY_TYPE_STRINGS[i].equals(displayType))
			{
				setDisplayType(i);
				return;
			}
		}
	}

    /*public boolean setDisplayVelocity(boolean display)
    {
        if(displayVelocity == display) return false;
        if(currentModel.win3d == null) return false;
        StsToggleButton toggleButton = (StsToggleButton)currentModel.win3d.getToolbarComponentNamed(StsVelocityAnalysisToolbar.NAME,  StsVelocityAnalysisToolbar.DISPLAYSTACK);
        if(toggleButton != null) toggleButton.setSelected(!display);
        displayVelocity = display;
        currentModel.viewObjectChanged(currentPreStackObject);
		currentModel.viewObjectRepaint(currentPreStackObject);
        return true;
    }

 public boolean getDisplayVelocity()
 {
	 return displayVelocity;
 }
*/
/*
    public void displayStack(boolean displayStack)
	{
		displayVelocity = !displayStack;
		currentModel.viewObjectRepaint(currentPreStackObject);
//        displayVelocity = true;
	}
*/
	/*
	 public void update()
	 {
	  ((StsPreStackLineSet)currentObject).update();
		}
	 */
/*
	public void setDisplayStack(boolean displayStack)
	{
        if(currentModel.win3d == null) return;
        StsToggleButton toggleButton = (StsToggleButton)currentModel.win3d.getToolbarComponentNamed(StsVelocityAnalysisToolbar.NAME,  StsVelocityAnalysisToolbar.DISPLAYSTACK);
        if(toggleButton != null) toggleButton.setSelected(displayStack);
        displayVelocity = !displayStack;
	}
*/
    public StsWiggleDisplayProperties getWiggleDisplayProperties()
    {
        if(wiggleDisplayProperties == null)
            wiggleDisplayProperties = new StsWiggleDisplayProperties(this, defaultWiggleDisplayProperties, "wiggleDisplayProperties");
        return wiggleDisplayProperties;
    }

    public StsWiggleDisplayProperties getPanelWiggleDisplayProperties()
    {
        if(panelWiggleDisplayProperties == null)
            panelWiggleDisplayProperties = new StsWiggleDisplayProperties(this, defaultWiggleDisplayProperties, "panelWiggleDisplayProperties");
        return panelWiggleDisplayProperties;
    }
}