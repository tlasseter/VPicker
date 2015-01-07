package com.Sts.DBTypes;

import com.Sts.Actions.Wizards.Color.*;
import com.Sts.DB.*;
import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.PreStack.*;
import com.Sts.Types.PreStack.StsSuperGather;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;
import com.Sts.Utilities.DataCube.*;
import com.Sts.Utilities.Interpolation.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.nio.*;
import java.text.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Nov 9, 2006
 * Time: 3:52:42 PM
 * To change this template use File | Settings | File Templates.
 */
abstract public class StsPreStackLineSet extends StsSeismicBoundingBox implements StsTreeObjectI, StsDialogFace, StsDataLineSetFace
{
    /** Predecessor to this volume */
	public StsPreStackLineSet previousVolume = null;
	/** velocity model generated in velocity analysis; used in NMOing volume */
	public StsPreStackVelocityModel velocityModel;
	/** input velocity model used in migration of this volume; used in DMOing volume */
	public StsPreStackVelocityModel inputVelocityModel;
	/** indicates lines are inlines; otherwise crosslines */
	public boolean isInline = true;
	/** current line being isVisible */
	public StsPreStackLine currentLine = null;
	/** array of lines in sorted order; taken from prestackLines initially read */
	public StsPreStackLine[] lines = null;
	/** All the colorscales being used */
	protected StsObjectRefList colorscales = null;
	/** The current colorscale viewed on the object panel */
	protected StsColorscale currentColorscale = null;
	/** indicates all lines/gathers in this data set are NMOed */
	public boolean isNMOed = false;
	/** display properties for all traces associated with this volume */
//    private StsWiggleDisplayProperties wiggleDisplayProperties = null;
	/** filter properties for all traces associated with this volume */
	public StsFilterProperties filterProperties = null;
	/** agc properties */
	public StsAGCPreStackProperties agcProperties = null;
	/** semblanceBytes display properties */
	public StsSemblanceDisplayProperties semblanceDisplayProperties = null;
	/** semblanceBytes edit properties */
	public StsSemblanceDisplayProperties semblanceEditProperties = null;
	/** semblanceBytes display properties */
	public StsSemblanceComputeProperties semblanceComputeProperties = null;
	/** velocity and z ranges for semblanceBytes */
	public StsSemblanceRangeProperties semblanceRangeProperties = null;
    /** datum properties */
    public StsDatumProperties datumProperties = null;
    /** superGather properties */
    public StsSuperGatherProperties superGatherProperties = null;
    /** CVS properties */
    public StsCVSProperties cvsProperties = null;

	/** range min of all trace offsets in volume */
	public float traceOffsetMin = StsParameters.largeFloat;
	/** range max of all trace offsets in volume */
	public float traceOffsetMax = -StsParameters.largeFloat;
	/** maximum number of traces in all gathers in volume */
	public int maxNTracesPerGather;
	/** minimum number of traces required to process analysis point */
//	int traceThreshold = 10;
	/** how many rows we skip between analysis points */
	int analysisRowInc = defaultAnalysisRowInc;
	/** how many cols we skip between analysis points */
	int analysisColInc = defaultAnalysisColInc;
    /** starting row number for analysis points */
    int analysisRowStart = 0;
    /** starting column number for analysis points */
    int analysisColStart = 0;
    /** clip percentage from medium for velocity profiles */
//    float corridorPercentage = defaultCorridorPercentage;
    /** Minumum number of traces in order to attempt a velocity picking */
    int traceThreshold = defaultTraceThreshold;

	/** semblanceBytes is in display or edit mode */
	public int semblanceMode = DISPLAY_MODE;
	/** maximum number of traces for all gathers in volume */
	public int nOffsetsMax;

	/** if true, display source positions on Z slice */
	boolean showSources = false;
	/** if true, display receiver positions on Z slice */
	boolean showReceivers = false;

    public String handVelName = null;
    public transient StsFile[] handVelFiles = null;

    /** Currently isVisible attribute */
    transient StsDisplaySeismicAttribute displayAttribute;
   /** Sources if they are available in attribute data */
    transient StsSeismicSources sources;
    /** Receivers if they are available in attribute data */
    transient StsSeismicReceivers receivers;
    /** Fold attribute data */
    transient StsPreStackFold fold;
    /** A fixed list of attributes which can be isVisible on Z slice */
	transient protected StsDisplaySeismicAttribute[] displayAttributes = new StsDisplaySeismicAttribute[0];

    transient StsProgressBarDialog progressBarDialog;
//    transient Runnable runExport;
//    transient Thread exportThread;

	transient public StsSuperGather[] superGathers = new StsSuperGather[0];
	transient public StsPreStackLineSetClass lineSetClass;
	transient StsColorList velocityColorList = null;
	/** colorscale for stacked data */
	transient protected StsColorscale seismicColorscale;
	/** colorscale for semblance data */
	transient protected StsColorscale semblanceColorscale;
   /** colorscale for semblance data */
	// transient protected StsColorscale stacksColorscale;

    transient VelocityProfileStatus velStat;
	transient protected boolean readoutEnabled = false;

	transient protected String displayCulture = DISPLAY_CULTURE_NONE;
	// the following are initialized after reading parameters file or reloading database
	transient protected StsSpectrumDialog spectrumDialog = null;
	transient protected boolean spectrumDisplayed = true;

//	transient protected int colorListNum = 0; // color list number for seismic data
//	transient protected boolean colorListChanged = true; // if colorscale has been edited, rebuild colorList displayList

//	transient boolean[] cursorDisplayChanged = new boolean[3];

	transient public StsColorList seismicColorList = null;
	transient public StsColorList semblanceColorList = null;
	// transient public StsColorList stacksColorList = null;

    /** pixels are isVisible as not-interpolated (isPixelMode == true) or as interpolated */
	transient boolean isPixelMode;
	/** indicates velocity is being isVisible; otherwise stack is isVisible */
//	transient boolean displayVelocity = true;

	transient byte stackOption;

    /** an iterator over the analysis grid; range is: analysisRowStart,analysisRowInc,analysisRowEnd  analysisColStart,analysisColInc,analysisColEnd */
    transient RowColIterator rowColIterator = null;

    static public final byte nullByte = StsParameters.nullByte;
    /** temporary mapped files used for storing float data along rows (lines) in this seismic data set */
    transient public StsCubeFileBlocks fileMapRowFloatBlocks;

	transient protected StsSeismicExportPanel exportPanel = null;

    static public byte[] byteTransparentTrace = null;
    static public float[] floatTransparentTrace = null;

    static StsPreStackLineSet[] lineSets = null;

	static private StsObjectPanel objectPanel = null;

	static public final int defaultAnalysisRowInc = 10;
	static public final int defaultAnalysisColInc = 10;
    static public final int defaultCorridorPercentage = 100;
    static public final int defaultTraceThreshold = 10;

    static public final byte CVS = 0;
    static public final byte SEMBLANCE = 1;

	static public final byte STACK_NONE = StsPreStackLineSet3dClass.STACK_NONE;
	static public final byte STACK_NEIGHBORS = StsPreStackLineSet3dClass.STACK_NEIGHBORS;
	static public final byte STACK_LINES = StsPreStackLineSet3dClass.STACK_LINES;
	static public final byte STACK_VOLUME = StsPreStackLineSet3dClass.STACK_VOLUME;

	static public final int DISPLAY_MODE = 0;
	static public final int EDIT_MODE = 1;
	static public final String DISPLAY_MODE_STRING = "Display Mode";
	static public final String EDIT_MODE_STRING = "Edit Mode";
	static public final String[] MODE_STRINGS = {DISPLAY_MODE_STRING, EDIT_MODE_STRING};

	static public String DISPLAY_CULTURE_NONE = "None";
	static public String[] displayCultures = new String[] {DISPLAY_CULTURE_NONE};

	static public final String SEISMIC_STRING = "Seismic";
	static public final String SEMBLANCE_STRING = "Semblance";
    static public final String STACKS_STRING = "Stacks";

    static public final int ATTRIBUTE_NONE = 0;
    static public final int ATTRIBUTE_FOLD = 1;
    static public final int ATTRIBUTE_SELEV = 2;
	static public final int ATTRIBUTE_RELEV = 3;
	static public final String ATTRIBUTE_NONE_STRING = "None";
	static public final String ATTRIBUTE_FOLD_STRING = "Fold";
	static public final String ATTRIBUTE_SELEV_STRING = "Source Elevation";
	static public final String ATTRIBUTE_RELEV_STRING = "Receiver Elevation";

    static public String NO_MODEL = "No Model";
    static public String HAND_VEL = "Hand Velocities";

	static protected StsComboBoxFieldBean colorscalesBean;
	static protected StsEditableColorscaleFieldBean colorscaleBean;
	static protected StsComboBoxFieldBean displayCultureBean;

	static public StsTimer timer = null;
	static public boolean runTimer = false;

	static public final boolean debug = false;
	static public final boolean debugTimer = false;

	abstract public boolean allocateVolumes(String mode, boolean loadFromFiles);
	abstract public int getVolumeRowColIndex(int row, int col);
	//  JKF 5/17/07 abstract public void setPlaneOK(boolean ok);
	abstract public void computePreStackVolume(byte displayType);
	abstract public StsPreStackLine getDataLine(int row, int col);
	abstract public Object[] getAvailableModelsList();
	abstract public StsPreStackVelocityModel constructVelocityModel();
	abstract public int getNColsForRow(int row);
	abstract public String getEmptyGatherDescription(int row, int col);
	abstract public String getFullGatherDescription(int row, int col, int nGatherTraces);
	abstract public void setInputVelocityModel(Object object);
	abstract public Object[] getAvailableModelsAndVolumesList();
	abstract public boolean isPlaneOK(byte displayType, int dir, int nPlane);
    abstract public void updatePreStackVolumes(StsRadialInterpolation interpolation);
    abstract public int[] adjustLimitRowCol(int nRow, int nCol);
    abstract public int[] getRowColFromCoors(float x, float y);

//	abstract public StsGather setCurrentDataRowCol(StsGLPanel3d glPanel3d, int row, int col);

	//abstract protected StsComboBoxFieldBean getColorscalesBean();
	//abstract protected StsEditableColorscaleFieldBean getColorscaleBean();
	//abstract protected StsComboBoxFieldBean getDisplayCultureBean();

    public StsPreStackLineSet()
	{
	}

	public StsPreStackLineSet(boolean persistent)
	{
		super(persistent);
    }

	static public void setStsObjectCopierCloneFlags(StsObjectCopier copier)
	{
		copier.setFieldCloneStsObject("previousVolume");
		copier.setFieldCloneStsObject("inputVelocityModel");
		copier.setFieldCloneStsObject("currentLine");
		copier.setFieldCloneStsObject("lines");
	}

	public boolean initialize(StsModel model)
	{
		try
		{
			lineSetClass = (StsPreStackLineSetClass) getStsClass();
//            lineSetClass.setCurrentObject(null);
//            if(velocityModel != null) velocityModel.setSeismicVolume(this);
//            if(inputVelocityModel != null) inputVelocityModel.setSeismicVolume(this);
//            setupColorscales();

			if (!allocateVolumes("rw", true)) return false;
            initializeColorscales();
            initializeDisplayAttributes();
            initCultures();
//			getWiggleDisplayProperties().buildBeanLists(this);
            initializePropertyPanels();
			setReceiver(showReceivers);
			setSource(showSources);

//			initializeDisplayAttribute();
			
			byte zDomainByte = StsParameters.getZDomainFromString(zDomain);
			zDomain = StsParameters.TD_ALL_STRINGS[zDomainByte];
			return model.getProject().checkSetZDomain(zDomainByte, zDomainByte);
		}
		catch (Exception e)
		{
			StsException.outputException("StsPreStackLineSet.classInitialize(StsModel) failed.", e, StsException.WARNING);
			return false;
		}
	}

    protected void checkSetOrigin(StsPreStackLine line)
    {
        if(originSet) return;
        checkSetOriginAndAngle(line.xOrigin, line.yOrigin, line.angle);
        isNMOed = line.isNMOed;
        isInline = line.isInline;
        byte zDomainByte = StsParameters.getZDomainFromString(line.zDomain);
        if (zDomainByte == StsParameters.TD_NONE)
            zDomain = StsParameters.TD_TIME_STRING;
        else
            zDomain = StsParameters.TD_ALL_STRINGS[zDomainByte];
        attributeNames = line.attributeNames;
        nAttributes = line.nAttributes;
        originSet = true;
    }

    private void initializePropertyPanels()
    {
        if(filterProperties != null) filterProperties.setParentObject(this);
        if(agcProperties != null) agcProperties.setParentObject(this);
        if(semblanceDisplayProperties != null) semblanceDisplayProperties.setParentObject(this);
        if(semblanceEditProperties != null) semblanceEditProperties.setParentObject(this);
        if(semblanceComputeProperties != null) semblanceComputeProperties.setParentObject(this);
        if(semblanceRangeProperties != null) semblanceRangeProperties.setParentObject(this);
        if(datumProperties != null) datumProperties.setParentObject(this);
        if(superGatherProperties != null) superGatherProperties.setParentObject(this);
        if(cvsProperties != null) cvsProperties.setParentObject(this);
        if(datumProperties != null) datumProperties.setParentObject(this);
    }

    private void initializeColorscales()
	{
        if(seismicColorscale == null)
            seismicColorscale = getColorscaleWithName(SEISMIC_STRING);
		if (seismicColorscale != null)
		{
			seismicColorList = new StsColorList(seismicColorscale);
			seismicColorscale.addActionListener(this);
		}
        if(semblanceColorscale == null)
            semblanceColorscale = getColorscaleWithName(SEMBLANCE_STRING);
		if (semblanceColorscale != null)
		{
			semblanceColorList = new StsColorList(semblanceColorscale);
			semblanceColorscale.addActionListener(this);
		}
    /*
        if(stacksColorscale == null)
            stacksColorscale = getColorscaleWithName(STACKS_STRING);
		if (stacksColorscale != null)
		{
			stacksColorList = new StsColorList(stacksColorscale);
			stacksColorscale.addActionListener(this);
		}
    */
    }

    public String getHandVelName() { return handVelName;}
    public void setHandVelName( String handVelName) { this.handVelName = handVelName;}
    public StsFile[] getHandVelFiles() { return handVelFiles;}
    public void setHandVelFiles( ArrayList<StsFile> files) 
    { 
        this.handVelFiles = new StsFile[files.size()];
        for (int i=0; i < files.size(); i++) handVelFiles[i] = files.get(i);
    }

    public StsPreStackLineSetClass getStsPreStackSeismicClass()
    {
        return (StsPreStackLineSetClass) getStsClass();
    }

	public StsColorscale getColorscaleWithName(String name)
	{
		int nColorscales = colorscales.getSize();
		for(int n = 0; n < nColorscales; n++)
		{
			StsColorscale colorscale = (StsColorscale)colorscales.getElement(n);
			if(colorscale.getName().equals(name)) return colorscale;
		}
		return null;
	}

    public StsColorscale[] getColorscaleList()
    {
        return (StsColorscale[])colorscales.getCastList(StsColorscale.class);
    }

    public void initializeSuperGathers()
    {
        superGathers = new StsSuperGather[0];
    }

   /** When db is reloaded, create these displayAttributes and initialize to current displayAttribute */
    protected void initializeDisplayAttributes()
    {
        StsNullAttribute none = new StsNullAttribute();
        addDisplayAttribute(none);
        fold = StsPreStackFold.constructor(this, currentModel);
        addDisplayAttribute(fold);
        sources = StsSeismicSources.constructor(this, currentModel);
        addDisplayAttribute(sources);
        receivers = StsSeismicReceivers.constructor(this, currentModel);
        addDisplayAttribute(receivers);
        displayAttribute = fold;
    }

    private void addDisplayAttribute(StsDisplaySeismicAttribute displayAttribute)
    {
        if(displayAttribute == null) return;
        displayAttributes = (StsDisplaySeismicAttribute[])StsMath.arrayAddElement(displayAttributes, displayAttribute, StsDisplaySeismicAttribute.class);
        addColorscale(displayAttribute.colorscale);
    }

	/*   public boolean classInitialize(StsModel model)
	  {
	 initializeStsClass();

	 if(!allocateVolumes("rw"))
	 {
	  return false;
	 }
	 initializeColorscale();
	 initCultures();
	 wiggleDisplayProperties.buildBeanLists(this);

	 setReceiver(showReceivers);
	 setSource(showSources);

//        String value = displayAttribute;
//        displayAttribute = DISPLAY_ATTRIBUTE_SEISMIC;
	 initializeDisplayAttribute();
	 byte zDomainByte = StsParameters.getZDomainFromString(zDomain);
	 zDomain = StsParameters.TD_ALL_STRINGS[zDomainByte];
	 if(!model.project.checkSetZDomain(zDomainByte, zDomainByte))
	 {
	  return false;
	 }
	 return true;
	   }*/

    protected void initialize()
    {
        this.lineSetClass = (StsPreStackLineSetClass) getStsClass();
        filterProperties = new StsFilterProperties(this, lineSetClass.defaultFilterProperties, "filterProperties");
        semblanceDisplayProperties = new StsSemblanceDisplayProperties(this, lineSetClass.defaultSemblanceDisplayProperties, StsSemblanceDisplayProperties.DISPLAY_MODE,
            "semblanceDisplayProperties");
        semblanceEditProperties = new StsSemblanceDisplayProperties(this, lineSetClass.defaultSemblanceEditProperties, StsSemblanceDisplayProperties.EDIT_MODE, "semblanceEditProperties");
        agcProperties = new StsAGCPreStackProperties(this, lineSetClass.defaultAGCProperties, "agcProperties");
        cvsProperties = new StsCVSProperties(this, lineSetClass.defaultCVSProperties, "cvsProperties");
        datumProperties = new StsDatumProperties(this, lineSetClass.defaultDatumProperties, "datumProperties");
        superGatherProperties = new StsSuperGatherProperties(this, lineSetClass.defaultSuperGatherProperties, "superGatherProperties");
        semblanceRangeProperties = new StsSemblanceRangeProperties(this, currentModel, "semblanceRangeProperties");
        semblanceComputeProperties = new StsSemblanceComputeProperties(this, lineSetClass.defaultSemblanceComputeProperties, "semblanceComputeProperties");
		setupColorscales();
//		setupDisplayAttributes();
		addToModel();
		semblanceRangeProperties.initializeSemblanceZRange(currentModel);
		currentModel.setCurrentObject(this);
		initialize(currentModel);
	}

    public StsPreStackVelocityModel getVelocityModel()
	{
		return velocityModel;
	}

    static public void checkAddToolbar(StsModel model, StsWin3dBase win3d)
    {
        if(model.classHasObjects(StsPreStackLineSet3d.class) || model.classHasObjects(StsPreStackLineSet2d.class))
            addToolbar(win3d);
    }

    static public void addToolbar()
	{
        addToolbar(currentModel.win3d);
	}

    static public void addToolbar(StsWin3dBase win3d)
    {
        if(win3d.hasToolbarNamed(StsSelectGatherToolbar.NAME)) return;
        StsSelectGatherToolbar toolbar = new StsSelectGatherToolbar(win3d);
        win3d.addToolbar(toolbar);
	}

    public String getGatherDescription(int row, int col, int nGatherTraces)
    {
        StsPreStackLine line = getDataLine(row, col);
        if (line == null || nGatherTraces == 0)
        {
            return getEmptyGatherDescription(row, col);
        }
        else
        {
            return getFullGatherDescription(row, col, nGatherTraces);
        }
    }

	public boolean isLineInSet(String stemname)
	{
		if (lines == null)
		{
			return false;
		}
		for (int j = 0; j < lines.length; j++)
		{
			if (lines[j] == null)
			{
				continue;
			}
			if (lines[j].stemname.equals(stemname))
			{
				return true;
			}
		}
		return false;
	}

	public int getNTracesInGather(int row, int col)
	{
		StsPreStackLine line = getDataLine(row, col);
		if (line == null)
		{
			return 0;
		}
		return line.getNTracesInGather(row, col);
	}

    public boolean hasGatherTraces(int row, int col)
    {
        return getNTracesInGather(row, col) > 0;
    }
    public boolean isNTracesBelowThreshold(int row, int col)
    {
        return getNTracesInGather(row, col) < traceThreshold;
    }

    /** Returns current gather set in the main window. The gather row&col corresponds to the main window cursor3d. */
    //TODO This is called by a several properties-related methods; which in turn call getNumberGathers() from this gather.
    //TODO Since there may be several window families with a different gather in each window (which presumably all must have the same number of gathers),
    //TODO it would make sense to get the number of gathers from the property itself directly rather than from a gather.
    public StsSuperGather getSuperGather()
	{
		return getSuperGather(currentModel.win3d);
	}

	public StsSuperGather getSuperGather(StsWin3dBase window)
	{
        if(window == null) return null;
        // if this is a parentWindow, it will be in the gatherList
		for (int n = 0; n < superGathers.length; n++)
		{
            StsWin3dBase gatherParentWindow = superGathers[n].parentWindow;
            if(gatherParentWindow == window || gatherParentWindow == window.parentWindow)
			{
				return superGathers[n];
			}
		}
		// need to create a superGather
		StsSuperGather superGather = StsSuperGather.constructor(currentModel, this, window.parentWindow);
        if(superGather == null) return null;
        // initialize it to zero-zero
//        superGather.initializeSuperGather(0, 0);
        superGathers = (StsSuperGather[]) StsMath.arrayAddElement(superGathers, superGather);
		return superGather;
	}

	public double[] getGatherAttribute(String attributeName, int row, int col)
	{
		try
		{
			StsPreStackLine line = getDataLine(row, col);
			return line.getAttributeArray(attributeName, row, col);
		}
		catch (Exception e)
		{
			StsException.outputException("StsPreStackLineSet.getGatherAttribute failed for attribute: " + attributeName + " row: " + row + " col: " + col, e, StsException.FATAL);
			return null;
		}
	}

    public void setIsNMOed(boolean isNMOed)
	{
        dbFieldChanged("isNMOed", isNMOed);
		this.isNMOed = isNMOed;
	}

	public boolean getIsNMOed()
	{
		return isNMOed;
	}

    public boolean flattenNMO()
    {
        return !isNMOed && lineSetClass.getFlatten();
    }

    public boolean canFlatten()
    {
        if(velocityModel == null) return false;
        return isNMOed || !isNMOed && velocityModel.hasProfiles();
    }

    public boolean canUnflatten()
    {
        if(velocityModel == null) return false;
        return  !isNMOed || isNMOed && velocityModel.hasProfiles();
    }

    public boolean hasVelocityProfiles()
    {
        return velocityModel != null && velocityModel.hasProfiles();
    }

    protected void setTraceOffsetRange(StsPreStackLine line)
	{
		maxNTracesPerGather = Math.max(maxNTracesPerGather, line.nOffsetsMax);
		traceOffsetMin = Math.min(traceOffsetMin, line.traceOffsetMin);
		traceOffsetMax = Math.max(traceOffsetMax, line.traceOffsetMax);
	}

	public boolean setGLColorList(GL gl, boolean nullsFilled, byte displayType, int shader)
	{
		StsColorList colorList = getSeismicColorList(displayType);
		if (colorList == null) return false;
		colorList.setGLColorList(gl, nullsFilled, shader);
		return true;
	}

	private StsColorList getSeismicColorList(byte displayType)
	{
        switch(displayType)
        {
            case StsPreStackLineSetClass.DISPLAY_VELOCITY:
                if(velocityModel != null)
                    return velocityModel.getVelocityColorList();
                else
                    return null;
            case StsPreStackLineSetClass.DISPLAY_STACKED:
			    return seismicColorList;
            case StsPreStackLineSetClass.DISPLAY_SEMBLANCE:
			    return semblanceColorList;
            default:
                return null;
        }
	}

	/** PostStack3d is not NMOed.  If we have interpreted profiles, use them to NMO (flatten) gather;
	 *  Otherwise if we have an inputVelocityModel (typically a segy volume has been read), then use it to NMO.
	 */
	public float[] getNMOVelocities(int row, int col)
	{
		if (velocityModel == null)
		{
			return null;
		}
		float[] velocities = velocityModel.getVelocities(row, col);
		return velocities;
	}

	public float[] getDNMOVelocities(int row, int col)
	{
		return null;
	}

	public float getTracesPerInch()
	{
		return getWiggleDisplayProperties().getTracesPerInch();
	}

	public float getInchesPerSecond()
	{
		return getWiggleDisplayProperties().getInchesPerSecond();
	}

    public byte getTimeUnits()
    {
        return currentModel.getProject().getTimeUnits();
    }

    public byte getDepthUnits()
    {
        return currentModel.getProject().getDepthUnits();
    }

    public byte getHorizontalUnits()
    {
        return currentModel.getProject().getXyUnits();
    }
	public boolean getReadoutEnabled()
	{
		return readoutEnabled;
	}

	public boolean getShowSource()
	{
		return showSources;
	}

	public boolean getShowReceiver()
	{
		return showReceivers;
	}

//    public boolean getShowVelStat() { return setShowVelStatusPoints; }

	public void setShowSource(boolean value)
	{
		if((currentLine.getAttributeIndex("SHT-X") == -1) || (currentLine.getAttributeIndex("SHT-Y") == -1))
		{
			new StsMessage(currentModel.win3d, StsMessage.WARNING, "Shot Positions are not available.");
			return;
		}
		if (setSource(value))
		{
			currentModel.win3dDisplayAll();
		}
	}

	protected boolean setSource(boolean value)
	{
		showSources = value;
		if (showSources && sources == null)
		{
			currentModel.addDisplayableInstance(sources);
			dbFieldChanged("showSources", showSources);
			currentModel.viewObjectRepaint(this, sources);
			return true;
		}
		else if (!showSources && sources != null)
		{
			currentModel.removeDisplayableInstance(sources);
			StsSeismicSourcesReceivers oldSources = sources;
			sources = null;
			dbFieldChanged("showSources", showSources);
			currentModel.viewObjectRepaint(this, oldSources);
			return true;
		}
		return false;
	}

	public void initializeVelStatusPoints()
	{
        boolean show = lineSetClass.getShowVelStat();
        setShowVelStatusPoints(show);
	}

    public void setShowVelStatusPoints(boolean show)
	{
		if (show && velStat == null)
		{
			velStat = new VelocityProfileStatus(this);
			currentModel.addDisplayableInstance(velStat);
		}
		else if (!show && velStat != null)
		{
			currentModel.removeDisplayableInstance(velStat);
			velStat = null;
		}
		return;
	}

	public void setShowReceiver(boolean value)
	{
		if ( (currentLine.getAttributeIndex(StsSEGYFormat.REC_X) == -1) || (currentLine.getAttributeIndex(StsSEGYFormat.REC_Y) == -1))
		{
			new StsMessage(currentModel.win3d, StsMessage.WARNING, "Receiver Position are not available.");
		}

		if (setReceiver(value))
		{
			currentModel.win3dDisplayAll();
		}

	}

	protected boolean setReceiver(boolean value)
	{
		showReceivers = value;
		if (showReceivers)
		{
			currentModel.addDisplayableInstance(receivers);
			dbFieldChanged("showReceivers", showReceivers);
			currentModel.viewObjectRepaint(this, receivers);
			return true;
		}
		else if (!showReceivers && receivers != null)
		{
			currentModel.removeDisplayableInstance(receivers);
			StsSeismicSourcesReceivers oldReceivers = receivers;
			receivers = null;
			dbFieldChanged("showReceivers", showReceivers);
			currentModel.viewObjectRepaint(this, oldReceivers);
			return false;
		}
		return false;
	}

	public void setReadoutEnabled(boolean enabled)
	{
		readoutEnabled = enabled;
		return;
	}

	public float getScaledValue(byte byteValue)
	{
		float f = (float) StsMath.signedByteToUnsignedInt(byteValue);
		return dataMin + (f / 254) * (dataMax - dataMin);
	}

	public void gatherPropertiesDialog()
	{
		gatherPropertiesDialog(currentModel.win3d);
	}

	public void gatherPropertiesDialog(Frame frame)
	{
		new StsOkApplyCancelDialog(frame, new StsDialogFace[]
								   {getWiggleDisplayProperties(), agcProperties, datumProperties, filterProperties}, "Gather Properties", false);
		// We will assume that some change has been made and all relevant views need to be redrawn.
		// The individual panelProperties will do whatever to prep for the redraw, but won't call a redraw.  So we may have unnecessary redraws,
		// but they shouldn't involve heavy recompute since textures will have been deleted and reconstructed just once.
		// viewObjectChangeRepaint() will only repaint those views who are affected by changes to this object
	}

	public void gatherDisplayPropertiesDialog()
	{
		gatherDisplayPropertiesDialog(currentModel.win3d);
	}

	public void gatherDisplayPropertiesDialog(Frame frame)
	{
		new StsOkApplyCancelDialog(frame, new StsDialogFace[]
								   {getWiggleDisplayProperties()}, "Gather Display Properties", false);
	}

	public void gatherAgcAndFilterPropertiesDialog()
	{
		gatherAgcAndFilterPropertiesDialog(currentModel.win3d);
	}

	public void gatherAgcAndFilterPropertiesDialog(Frame frame)
	{
		new StsOkApplyCancelDialog(frame, new StsDialogFace[]
								   {agcProperties, filterProperties}, "Gather AGC & Filter Properties", false);
	}

	public void semblancePropertiesDialog(Frame frame, StsMouse mouse)
	{
		StsBatchOkApplyCancelDialog dialog = new StsBatchOkApplyCancelDialog(frame, this, this.SEMBLANCE, new StsDialogFace[]
								   {semblanceComputeProperties, semblanceRangeProperties, getSemblanceDisplayProperties(),
                                    superGatherProperties, agcProperties, datumProperties, filterProperties},
								   "Semblance Properties", false);
		dialog.setLocation(mouse.getX(), mouse.getY());
	}

	public void semblanceDisplayPropertiesDialog()
	{
		semblanceDisplayPropertiesDialog(currentModel.win3d);
	}

	public void semblanceDisplayPropertiesDialog(Frame frame)
	{
		new StsOkApplyCancelDialog(frame, new StsDialogFace[]
								   {semblanceComputeProperties, semblanceRangeProperties, getSemblanceDisplayProperties()}, "Semblance Display Properties", false);
	}

	public void semblanceAgcAndFilterPropertiesDialog()
	{
		semblanceAgcAndFilterPropertiesDialog(currentModel.win3d);
	}

	public void semblanceAgcAndFilterPropertiesDialog(Frame frame)
	{
		new StsOkApplyCancelDialog(frame, new StsDialogFace[]
								   {agcProperties, filterProperties}, "Semblance AGC and Filter Properties", false);
	}

    public void cvsDisplayPropertiesDialog()
    {
        cvsPropertiesDialog(currentModel.win3d, null);
	}

	public void cvsPropertiesDialog(Frame frame, StsMouse mouse)
	{
		StsBatchOkApplyCancelDialog dialog = new StsBatchOkApplyCancelDialog(frame, this, this.CVS, new StsDialogFace[]
		   {cvsProperties, semblanceRangeProperties, getSemblanceDisplayProperties(), getWiggleDisplayProperties(), superGatherProperties,
                   agcProperties, filterProperties}, "CVS/VVS Properties", false);
		if (mouse != null) 
		{
		    dialog.setLocation(mouse.getX(), mouse.getY());
		}
    }

	public StsSemblanceDisplayProperties getSemblanceDisplayProperties()
	{
		if (semblanceMode == DISPLAY_MODE)
		{
			return semblanceDisplayProperties;
		}
		else
		{
			return semblanceEditProperties;
		}
	}

    public void displayWiggleProperties()
    {
        new StsOkApplyCancelDialog(currentModel.win3d, new StsDialogFace[] {getWiggleDisplayProperties()}, "Wiggle Properties", false);
    }
    public void displaySemblanceProperties()
    {
        new StsOkApplyCancelDialog(currentModel.win3d, new StsDialogFace[] {getSemblanceDisplayProperties() }, "Velocity Display Properties", false);
    }
    public void displayAGCProperties()
    {
        new StsOkApplyCancelDialog(currentModel.win3d, new StsDialogFace[] {agcProperties}, "AGC Properties", false);
    }
    public void displayFilterProperties()
    {
        new StsOkApplyCancelDialog(currentModel.win3d, new StsDialogFace[] {filterProperties}, "Filter Properties", false);
    }
    public void displayCVSProperties()
    {
        new StsOkApplyCancelDialog(currentModel.win3d, new StsDialogFace[] {cvsProperties}, "Edit CVS/VVS Properties", false);
    }
    public void displayDatumProperties()
    {
        new StsOkApplyCancelDialog(currentModel.win3d, new StsDialogFace[] {datumProperties}, "Edit Datum Properties", false);
    }
    public void displaySuperGatherProperties()
    {
        new StsOkApplyCancelDialog(currentModel.win3d, new StsDialogFace[] {superGatherProperties}, "Edit Super Gather Properties", false);
    }
	public StsSemblanceComputeProperties getSemblanceComputeProperties()
	{
		return semblanceComputeProperties;
	}
    
    public double getDatumShift()
    {
        double datumShift = 0.0;
        if( datumProperties != null)
            datumShift = datumProperties.getDatum();
        return datumShift;
    }

	public boolean agcOrFilterPropertiesChanged(Object object)
	{
        if(object instanceof StsAGCPreStackProperties) return true;
        if(object instanceof StsFilterProperties) return true;
        if(object instanceof StsDatumProperties) return true;
        return false;
	}
/*
	public boolean cvsPropertiesChanged()
	{
		return cvsProperties.isChanged() || semblanceRangeProperties.isChanged() || agcProperties.isChanged() ||
											superGatherProperties.isChanged() || filterProperties.isChanged();
	}
*/
	protected void checkMemoryStatus(String message)
	{
		NumberFormat numberFormat = NumberFormat.getInstance();
		String freeMemory = numberFormat.format(Runtime.getRuntime().freeMemory());
		String totalMemory = numberFormat.format(Runtime.getRuntime().totalMemory());
		String maxMemory = numberFormat.format(Runtime.getRuntime().maxMemory());
		if (debug)
		{
			System.out.println(message + ". freeMemory: " + freeMemory +
							   " totalMemory: " + totalMemory + " maxMemory: " + maxMemory);
		}
	}

	protected void setupColorscales()
	{
		try
		{
			if (seismicColorscale == null)
			{
				StsSpectrumClass spectrumClass = currentModel.getSpectrumClass();
				StsSpectrum seismicSpectrum = spectrumClass.getSpectrum(lineSetClass.getSeismicSpectrumName());
				seismicColorscale = new StsColorscale(SEISMIC_STRING, seismicSpectrum, -1.0f, 1.0f);
			}
			addColorscale(seismicColorscale);

			if (semblanceColorscale == null)
			{
				StsSpectrumClass spectrumClass = currentModel.getSpectrumClass();
				StsSpectrum semblanceSpectrum = spectrumClass.getSpectrum(lineSetClass.getSemblanceSpectrumName());
				semblanceColorscale = new StsColorscale(SEMBLANCE_STRING, semblanceSpectrum, 0.0f, 0.25f);
                semblanceColorscale.setEditMax(0.25f);
                semblanceColorscale.setCompressionMode(StsColorscale.COMPRESSED);
                semblanceColorscale.setTransparencyMode(false);
			}
			addColorscale(semblanceColorscale);
        /*
            if (stacksColorscale == null)
			{
				StsSpectrumClass spectrumClass = currentModel.getSpectrumClass();
				StsSpectrum stacksSpectrum = spectrumClass.getSpectrum(lineSetClass.getStacksSpectrumName());
				stacksColorscale = new StsColorscale(STACKS_STRING, stacksSpectrum, -1.0f, 1.0f);
			}
			addColorscale(stacksColorscale);
	     */
        }
		catch (Exception e)
		{
			StsException.outputException("StsPreStackLineSet.initializeColorscale() failed.", e, StsException.WARNING);
		}
	}

	public void addColorscale(StsColorscale colorscale)
	{
		if(colorscales == null) colorscales = StsObjectRefList.constructor(2, 2, "colorscales", this);
		colorscales.add(colorscale);
//		colorscales = (StsColorscale[]) StsMath.arrayAddElementNoRepeat(colorscales, colorscale);
		setCurrentColorscale(colorscale);
		// if(colorscalesBean == null) return;
		// colorscalesBean.setListItems(colorscales.getElements());
		// colorscalesBean.setSelectedItem(colorscale);
	}

	public void setCurrentColorscale(StsColorscale colorscale)
	{
		if (currentColorscale == colorscale)
		{
			return;
		}
		currentColorscale = colorscale;
		if(isPersistent()) dbFieldChanged("currentColorscale", currentColorscale);
		if(colorscaleBean != null) colorscaleBean.setValueObject(colorscale);
		if(colorscalesBean != null) colorscalesBean.setSelectedItem(colorscale);
	}

	public StsColorscale getCurrentColorscale()
	{
		return currentColorscale;
	}

	public StsColorscale getSeismicColorscale()
	{
//		setDataHistogram();
		return seismicColorscale;
	}
/*
    public StsColorscale getStacksColorscale()
	{
		return stacksColorscale;
	}
*/
	/*
	 public void setDataHistogram()
	 {
	  getColorscaleBean().setHistogram(dataHist);
	 }*/

	public void setSeismicColorscale(StsColorscale colorscale)
	{
		this.seismicColorscale = colorscale;
		currentModel.viewObjectRepaint(this, colorscale);
	}
/*
	public void setStacksColorscale(StsColorscale colorscale)
	{
		this.stacksColorscale = colorscale;
		currentModel.viewObjectRepaint(this, colorscale);
	}
*/
    public void clearAllGathers()
	{
		for (int n = 0; n < lines.length; n++)
		{
			if (lines[n] != null)
			{
				lines[n].clearFileBlocks();
			}
		}
	}
/*
    public int getAnalysisRowNumStart()
    {
        return getRowNumFromRow(analysisRowStart);
    }

    public int getAnalysisColNumStart()
    {
        return getColNumFromCol(analysisColStart);
    }
*/

    public float getAnalysisColNumStart()
    {
        return getColNumFromCol(analysisColStart);
    }

    public void setAnalysisRowNumStart(float rowStartNum)
    {
        analysisRowStart = getRowFromRowNum(rowStartNum);
        dbFieldChanged("analysisRowStart", analysisRowStart);
        currentModel.viewObjectRepaint(this, this);
    }

    public float getAnalysisRowNumStart()
    {
        return getRowNumFromRow(analysisRowStart);
    }

    public void setAnalysisColNumStart(float colStartNum)
    {
        analysisColStart = getColFromColNum(colStartNum);
        dbFieldChanged("analysisColStart", analysisColStart);
        currentModel.viewObjectRepaint(this, this);
    }

    public float getAnalysisRowNumInc()
    {
        return analysisRowInc*rowNumInc;
    }

    public void setAnalysisRowNumInc(float analysisRowNumInc)
    {
        analysisRowInc = Math.round(analysisRowNumInc/rowNumInc);
        dbFieldChanged("analysisRowInc", analysisRowInc);
        currentModel.viewObjectRepaint(this, this);
    }

    public float getAnalysisColNumInc()
    {
        return analysisColInc*colNumInc;
    }

    public void setAnalysisColNumInc(float analysisColNumInc)
    {
        analysisColInc = Math.round(analysisColNumInc/colNumInc);
        dbFieldChanged("analysisColInc", analysisColInc);
        currentModel.viewObjectRepaint(this, this);
    }


    public void setAnalysisRowStart(int rowStart)
    {
        analysisRowStart = rowStart;
        dbFieldChanged("analysisRowStart", analysisRowStart);
        currentModel.viewObjectRepaint(this, this);
    }

    public int getAnalysisRowStart()
    {
        return analysisRowStart;
    }

    public void setAnalysisColStart(int analysisColStart)
    {
        this.analysisColStart = analysisColStart;
        dbFieldChanged("analysisColStart", analysisColStart);
        currentModel.viewObjectRepaint(this, this);
    }

    public int getAnalysisColStart()
    {
        return analysisColStart;
    }

    public int getAnalysisColInc()
    {
        return analysisColInc;
    }
    public int getAnalysisRowInc()
    {
        return analysisRowInc;
    }
    public void setAnalysisColInc(int analysisColInc)
    {
        this.analysisColInc = analysisColInc;
        dbFieldChanged("analysisColInc", analysisColInc);
        currentModel.viewObjectRepaint(this, this);
    }
    public void setAnalysisRowInc(int analysisRowInc)
    {
        this.analysisRowInc = analysisRowInc;
        dbFieldChanged("analysisRowInc", analysisRowInc);
        currentModel.viewObjectRepaint(this, this);
    }
    public void actionPerformed(ActionEvent e)
	{
		if (! (e.getSource() instanceof StsColorscale))
			return;
//		int eventID = e.getID();
//        if(eventID == MouseEvent.MOUSE_DRAGGED) return;
		StsColorscale colorscale = (StsColorscale) e.getSource();
		if (semblanceColorscale == colorscale)
		{
			semblanceColorList.setColorListChanged(true);
            currentModel.clearDisplayTextured3dCursors(this);
			currentModel.viewObjectRepaint(colorscale, semblanceColorList);
		}
		else if(seismicColorscale == colorscale)
		{
			seismicColorList.setColorListChanged(true);
            currentModel.clearDisplayTextured3dCursors(this);
            currentModel.viewObjectRepaint(colorscale, seismicColorList);
		}
        /*
        else if(stacksColorscale == colorscale)
		{
			stacksColorList.setColorListChanged(true);
//			currentModel.viewObjectChanged(semblanceColorList);
			currentModel.viewObjectChangedAndRepaint(this, stacksColorList);
		}
		*/
    }

	public void semblanceColorListChanged()
	{
		semblanceColorList.setColorListChanged(true);
		currentModel.viewObjectRepaint(this, semblanceColorList);
	}

	public void seismicColorListChanged()
	{
		seismicColorList.setColorListChanged(true);
		currentModel.viewObjectRepaint(this, seismicColorList);
	}
/*
	public void stacksColorListChanged()
	{
		stacksColorList.setColorListChanged(true);
		currentModel.viewObjectRepaint(this, stacksColorList);
	}
*/
    /** volume has been removed from display; delete allocated planes, etc */
	public void deleteTransients()
	{
		clearCache();
	}

	public void clearCache()
	{
		StsBlocksMemoryManager memoryManager = currentModel.getProject().getBlocksMemoryManager();
		memoryManager.clearAllBlocks();
	}

	public boolean stackOptionChanged()
	{
		byte currentStackOption = lineSetClass.getStackOption();
		if (currentStackOption == stackOption)
		{
			return false;
		}
		stackOption = currentStackOption;
		return true;
	}

	public void setStackOption(byte option)
	{
		stackOption = option;
		// JKF setPlaneOK(false);
	}

	public boolean isPixelModeChanged()
	{
		if (lineSetClass.isPixelMode == isPixelMode)
		{
			return false;
		}
		isPixelMode = !isPixelMode;
		return true;
	}

	public void setMode(int mode)
	{
		this.semblanceMode = mode;
	}

	public void setDisplayCulture(String attribute)
	{
		if (displayCulture == attribute)
		{
			return;
		}
		displayCulture = attribute;
//		cursorDisplayChanged[ZDIR] = true;
		currentModel.win3dDisplayAll();
	}

	public String getDisplayCulture()
	{
		return displayCulture;
	}

	public void initCultures()
	{
		boolean found = false;
		if (lines == null)
		{
			return; // can happen if 2d only
		}
		for (int k = 0; k < lines.length; k++)
		{
			if (lines[k] == null)
			{
				continue;
			}
			String[] names = lines[k].getAttributes();
			for (int i = 0; i < names.length; i++)
			{
				found = false;
				for (int j = 0; j < displayCultures.length; j++)
				{
					if (names[i].equals(displayCultures[j]))
					{
						found = true;
						break;
					}
				}
				if (!found)
				{
					displayCultures = (String[]) StsMath.arrayAddElement(displayCultures, names[i]);
				}
			}
		}
//        getDisplayCultureBean().setListItems(displayCultures);
	}

	public String getSegyFileDate()
	{
		if (segyLastModified == 0)
		{
			File segyFile = new File(segyDirectory + segyFilename);
			if (segyFile != null)
			{
				segyLastModified = segyFile.lastModified();
			}
		}
		DateFormat dateFormat = DateFormat.getDateTimeInstance();
		return dateFormat.format(new Date(segyLastModified));
	}

    public StsDisplaySeismicAttribute[] getDisplayAttributes() { return displayAttributes; }

    public StsDisplaySeismicAttribute getDisplayAttribute()
	{
		return displayAttribute;
	}
    //TODO save the string as a object.field property
	public void setDisplayAttribute(StsDisplaySeismicAttribute displayAttribute)
	{
        if(displayAttributes == null) return;
        if(displayAttribute == this.displayAttribute) return;
        this.displayAttribute = displayAttribute;
        displayAttribute.colorListChanged();
        currentModel.viewObjectChangedAndRepaint(this, this);
    }

    public boolean delete()
	{
		close();
		// the class could be either StsSeismicVolume or StsVirtualVolume
		StsClass stsClass = currentModel.getStsClass(getClass());
		StsObject currentObject = stsClass.getCurrentObject();
		stsClass.delete(this);
		super.delete(); // sets this.index to -1

		if (currentObject == this)
		{
			int nVolumes = stsClass.getSize();
			if (nVolumes > 0)
			{
				currentObject = stsClass.getElement(nVolumes - 1);
			}
			else
			{
				currentObject = null;

			}
			stsClass.setCurrentObject(currentObject);
            currentModel.toggleOnCursor3dObject(currentObject);
		}
		return true;
	}

	public int getTraceThreshold()
	{
        //if(velocityModel == null) return defaultTraceThreshold;   //this is causing many problems, when defining velocity analysis parameters velocityModel will be null!!!
		return traceThreshold;
	}

	public void setTraceThreshold(int threshold)
	{
        //if(velocityModel == null) return;
        if(threshold == traceThreshold) return;
        traceThreshold = threshold;
        fieldChanged("traceThreshold", traceThreshold);
        currentModel.repaintViews(StsView3d.class);
    }
/*
    public float getCorridorPercentage()
    {
        if(velocityModel == null) return defaultCorridorPercentage;
        return corridorPercentage;
    }

    public void setCorridorPercentage(float percentage)
    {
        if(velocityModel == null) return;
        corridorPercentage = percentage;
    }
*/
    public void previousProfile()
    {
        previousProfile(currentModel.win3d);
    }

	public void previousProfile(StsWin3dBase window)
	{
        moveToProfile(window, -1);
    }

    public void nextProfile()
    {
        nextProfile(currentModel.win3d);
    }
	public void nextProfile(StsWin3dBase window)
	{
        moveToProfile(window, 1);
    }

    public void moveToProfile(StsWin3dBase window, int direction)
	{
//		updateCurrentVelocityProfile();
        if(window == null) window = currentModel.win3d;
        StsCursor3d cursor = window.getCursor3d();
        if (cursor == null) return;
        int row = cursor.getRow();
        int col = cursor.getCol();
        
        byte profileOption = lineSetClass.getProfileOption();
		int[] rowCol;
		int nProfiles = 0;
		if(velocityModel != null)
		{
			nProfiles = velocityModel.velocityProfiles.getSize();
		}
		if (profileOption == StsPreStackLineSet3dClass.EVERY_PROFILE || profileOption == StsPreStackLineSet3dClass.UNPICKED_PROFILES && nProfiles == 0)
		{
            // this initializes iterator to current position and sets next to this row-col
            // so advance it, so iterator will actually return the next one on the next() call
            rowColIterator = getRowColIterator(row, col, direction);
            if(rowColIterator.hasNext()) rowColIterator.next();
            if (debug) System.out.println("Row=" + row + " Col=" + col);

            int tempTraceThreshold = traceThreshold;
            while(tempTraceThreshold >= 1)
            {
                while(rowColIterator.hasNext())
                {
                    rowCol = (int[])rowColIterator.next();
                    if(getNTracesInGather(rowCol[0], rowCol[1]) >= tempTraceThreshold)
                    {
                        jumpToRowCol(rowCol, window);
                        return;
                    }
                }
                tempTraceThreshold--;
            }
		}
		else if (profileOption == StsPreStackLineSet3dClass.PICKED_PROFILES)
		{
            if(direction > 0)
                rowCol = getNextPickedProfileRowCol(row, col);
            else
                rowCol = getPrevPickedProfileRowCol(row, col);
			jumpToRowCol(rowCol, window);
		}
		else // UNPICKED_PROFILES
		{
           // this initializes iterator to current position and sets next to this row-col
            // so advance it, so iterator will actually return the next one on the next() call
            rowColIterator = getRowColIterator(row, col, direction);
            if(rowColIterator.hasNext()) rowColIterator.next();

            int tempTraceThreshold = traceThreshold;
            while(tempTraceThreshold >= 1)
            {
                while(rowColIterator.hasNext())
                {
                    rowCol = (int[])rowColIterator.next();
                    if(!hasProfile(rowCol) && getNTracesInGather(rowCol[0], rowCol[1]) >= tempTraceThreshold)
                    {
                        jumpToRowCol(rowCol, window);
                        return;
                    }
                }
                tempTraceThreshold--;
                rowColIterator.initialize(row, col, direction);
            }
        }
	}

	public void jumpToRowCol(int[] rowCol, StsWin3dBase window)
	{
		if (rowCol == null) return;
		int row = rowCol[0];
		int col = rowCol[1];
        StsSuperGather superGather = setCurrentDataRowCol(window, row, col);
        float[] xy = getXYCoors(row, col);
        currentModel.win3d.adjustSlider(YDIR, xy[1]);
        currentModel.win3d.adjustSlider(XDIR, xy[0] ); //these methods adjust both slider and cursor (indirectly via "stateChanged")
        currentModel.win3d.adjustCursorXY(xy[0], xy[1]); //for 2d, slider is inaccurate (off by up to "Step" distance), so need to reset cursor to true coordinates
        currentModel.win3d.cursor3d.setRow(row);
        currentModel.win3d.cursor3d.setCol(col);
        currentModel.viewObjectChangedAndRepaint(this, superGather);
    }

    /** row or col may be -1 meaning this index should be left the same as currently */
	public StsSuperGather setCurrentDataRowCol(StsWin3dBase window, int row, int col)
	{
	    if(velocityModel != null) 
	    {
	        velocityModel.checkCurrentVelocityProfile(window);
	        // float[] xy = getXYCoors(row, col);
	        // velocityModel.setCurrentXY(xy[0], xy[1]);
	    }
        StsSuperGather gather = getSuperGather(window);
        if(gather == null)  return null;
        // gather.classInitialize(...) returns true if already initialized
        if (!gather.initializeSuperGather(row, col)) return null;
		if (lines == null) return null;
		currentLine = getDataLine(gather.superGatherRow, gather.superGatherCol);
        if(debug)
        {
            if(currentLine == null)
                System.out.println("StsVelocityAnalysisEdit.line2dTextureChanged(): line2d is null");
            else
                System.out.println("StsVelocityAnalysisEdit.line2dTextureChanged(): line2d is " + currentLine.lineIndex);
        }
        if (currentLine == null) return null;
        // If stack or semblance is being isVisible, switch to velocity display unless we are only displaying NEIGHBORHOOD.
        if(velocityModel != null)
        {
            byte displayType =  lineSetClass.getDisplayType();
            if(displayType != StsPreStackLineSetClass.DISPLAY_VELOCITY)
            {
                boolean updateData = !isPlaneOK(displayType, StsRotatedGridBoundingBox.XDIR, col) || !isPlaneOK(displayType, StsRotatedGridBoundingBox.YDIR, row);
                if(updateData && lineSetClass.getStackOption() == StsPreStackLineSetClass.STACK_LINES)
                   lineSetClass.setDisplayType(StsPreStackLineSetClass.DISPLAY_VELOCITY);
            }
        }
        return gather;
	}

    public Iterator getGatherIterator()
    {
        return new StsTestGatherIterator(this);
    }

    /** A pre-stack seismic volume may have been created from the migration of a previous volume so we
     * need to know the predecessor volume.
     */
    public StsPreStackLineSet getPredecessor() { return previousVolume; }
    public void setPredecessor(StsPreStackLineSet volume)
    {
        previousVolume = volume;
    }
        /** given a row-col which may or may not be on a regular grid (rows: 0, 0+rowInc.. cols: 0, 0+colInc...),
	 *  go to a next (direction == 1) or prev (direction == 0) row-col which is on the regular grid.
	 *  First round row-col to lower-left row-col on grid.  Compute the index of this point on the regular grid.
	 *  For a row-ordered grid (isInline==true), index = row*nCols + nCol.
	 *  If we want next row-col, increment index by 1; if we want prev row-col and original point was not on grid,
	 *  decrement index by 1.  If we wanted prev and index below 0, set index to max.
	 *  If we wanted next and index above max, set index to 0.
	 *  Now convert from grid row-col back to normal row-col.
	 */

	public int[] getPrevRowCol(int row, int col)
	{
		return advanceRowCol(row, col, -1);
	}

	public int[] getNextRowCol(int row, int col)
	{
		return advanceRowCol(row, col, 1);
	}


    public int[] advanceRowCol(int row, int col, int direction)
    {
        rowColIterator = getRowColIterator(row, col, direction);
        return (int[])rowColIterator.next();
    }

    public int getRowColIndex(int[] rowCol)
	{
		return getVolumeRowColIndex(rowCol[0], rowCol[1]);
	}

    // Returns the row, col of the next profile which may be the row, col past in.
    private boolean hasProfile(int[] rowCol)
    {
        StsObjectRefList velocityProfiles = velocityModel.velocityProfiles;
        int nProfiles = velocityProfiles.getSize();
        if (nProfiles == 0)
        {
            return false;
        }
        int currentRowColIndex = getVolumeRowColIndex(rowCol[0], rowCol[1]);
        for (int n = 0; n < velocityProfiles.getSize(); n++)
        {
            StsVelocityProfile profile = (StsVelocityProfile) velocityProfiles.getElement(n);
            int profileRowColIndex = getVolumeRowColIndex(profile.row, profile.col);
            if(currentRowColIndex == profileRowColIndex)
                return true;
        }
        return false;
    }
	private int[] getNextPickedProfileRowCol(int row, int col)
	{
		StsObjectRefList velocityProfiles = velocityModel.velocityProfiles;
        if(velocityProfiles == null) return null;
		int nProfiles = velocityProfiles.getSize();
		if (nProfiles == 0)
		{
			return null;
		}
		int currentRowColIndex = getVolumeRowColIndex(row, col);
		int nextRowColIndex = -1;
		for (int n = 0; n < velocityProfiles.getSize(); n++)
		{
			StsVelocityProfile profile = (StsVelocityProfile) velocityProfiles.getElement(n);
			int prevRowColIndex = nextRowColIndex;
			nextRowColIndex = getVolumeRowColIndex(profile.row, profile.col);
			if (currentRowColIndex >= prevRowColIndex && currentRowColIndex < nextRowColIndex)
			{
				return new int[]
					{profile.row, profile.col};
			}
		}
		// must be > last, so return the first
		StsVelocityProfile profile = (StsVelocityProfile) velocityProfiles.getElement(0);
		return new int[]
			{profile.row, profile.col};
	}

	private int[] getPrevPickedProfileRowCol(int row, int col)
	{
		StsObjectRefList velocityProfiles = velocityModel.velocityProfiles;
		int nProfiles = velocityProfiles.getSize();
		if (nProfiles == 0)
		{
			return null;
		}
		int currentRowColIndex = getVolumeRowColIndex(row, col);
		int prevRowColIndex = Integer.MAX_VALUE;
		for (int n = nProfiles - 1; n >= 0; n--)
		{
			StsVelocityProfile profile = (StsVelocityProfile) velocityProfiles.getElement(n);
			int nextRowColIndex = getVolumeRowColIndex(profile.row, profile.col);
			if (currentRowColIndex > nextRowColIndex && currentRowColIndex <= prevRowColIndex)
			{
				return new int[] {profile.row, profile.col};
			}
		}
		// must be < first, so return the last
		StsVelocityProfile profile = (StsVelocityProfile) velocityProfiles.getElement(nProfiles - 1);
		return new int[] {profile.row, profile.col};
	}

	 public void objectPropertiesChanged(Object object)
	 {
		if(!agcOrFilterPropertiesChanged(object)) return;
        currentModel.viewObjectChangedAndRepaint(this, this);
     }

	public void objectChanged()
	{
	}

    public boolean getApplyFilter(byte dataType)
    {
        if(filterProperties.getApplyBoxFilter(dataType)) return true;
        if(filterProperties.getApplyBWFilter(dataType)) return true;
		if(agcProperties.getApplyAGC(dataType)) return true;
        return false;
    }

    public void setIsVisible(boolean isVisible)
	{
		super.setIsVisible(isVisible);
        currentModel.viewObjectRepaint(this, this);
    }

    public void setVelocityModel(StsPreStackVelocityModel newVelocityModel)
	{
		if (newVelocityModel == null) return;
		dbFieldChanged("velocityModel", velocityModel, newVelocityModel);
        velocityModel = newVelocityModel;
    }

	public Object getInputVelocityModel()
	{
		if (inputVelocityModel == null)
		{
			return NO_MODEL;
		}
		else
		{
			return inputVelocityModel;
		}
	}

	public boolean hasExportableData(int row)
	{
        if(velocityModel == null) return false;
        return velocityModel.interpolation.rowHasDataPoint(row);
	}

	public String getCurrentLineFilename()
	{
		if(currentLine != null)return currentLine.stemname;
		return "";
	}

	public Object[] getChildren()
	{
		return new Object[0];
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
		lineSetClass.selected(this);
		currentModel.getGlPanel3d().checkAddView(StsView3d.class);
//        currentModel.glPanel3d.cursor3d.setCurrentSeismicVolume(this);
        // colorscalesBean.setListItems(colorscales.getElements());
		currentColorscale = (StsColorscale)colorscales.getElement(0);
		currentModel.win3dDisplayAll();
	}

	public boolean anyDependencies()
	{
		return false;
	}


    static public void checkTransparentTrace(int nSlices)
    {
        if(byteTransparentTrace != null && byteTransparentTrace.length > nSlices) return;
        byteTransparentTrace = new byte[nSlices];
        floatTransparentTrace = new float[nSlices];
        for(int n = 0; n < nSlices; n++)
        {
            byteTransparentTrace[n] = -1;
            floatTransparentTrace[n] = 0;
        }
        // TODO: classInitialize floatTransparentTrace to appropriate value
    }

    public StsWiggleDisplayProperties getWiggleDisplayProperties()
    {
        return lineSetClass.getWiggleDisplayProperties();
    }

    public StsCVSProperties getCVSProperties()
    {
        return cvsProperties;
    }

    public String getDescription()
    {
 		return new String("Name: " + stemname + "\n" + "SegY File: " + segyFilename + "\n" +
            "SegY Last Modified On: " + new Date(segyLastModified).toString() + "\n" +
            "Data Range: " + String.valueOf(dataMin) + " to " + String.valueOf(dataMax) + "\n" +
            "X Origin: " + String.valueOf(xOrigin) + " Y Origin: " + String.valueOf(yOrigin));
    }

    public StsSpiralRadialInterpolation getSpiralRadialInterpolation(int nMinPoints, int nMaxPoints)
    {
        return new StsSpiralRadialInterpolation(this, nMinPoints, nMaxPoints);
    }

	public ByteBuffer getAttributeByteBuffer(int dir)
	{
		if(dir != ZDIR) return null;
        if(displayAttribute.byteBuffer == null)
            displayAttribute.computeBytes();
        return displayAttribute.byteBuffer;
	}

    public boolean hasAttributeData()
    {
        return displayAttribute.getName() != ATTRIBUTE_NONE_STRING;
    }

    public boolean isAttributeTextureChanged()
    {
        return displayAttribute.byteBuffer == null;
    }

    class VelocityProfileStatus implements StsInstance3dDisplayable
	{
		StsPreStackVelocityModel velocityVolume = null;
		StsPreStackLineSet volume = null;

//		StsColor HAND_VELOCITY = StsColor.GREEN;
//		StsColor PLANNED_HAND_VELOCITY = StsColor.YELLOW;

		public VelocityProfileStatus(StsPreStackLineSet vol)
		{
			velocityVolume = vol.velocityModel;
			volume = vol;
		}

		public void display(StsGLPanel3d glPanel3d)
		{
			if (velocityModel == null) return;
		    if(!velocityModel.getDisplayAnalysisPoints()) return;
            
            int size = lineSetClass.getDisplayPointSize(); //currentModel.getIntProperty("displayPointSize");
            int buff = 4;
			try
			{
			    float[] xyz = currentModel.win3d.getCursor3d().getCurrentCoordinates();
	            //xyz[2] = currentModel.win3d.getCursor3d().getCurrentDirCoordinate(StsCursor3d.ZDIR);
                StsSuperGather gather = getSuperGather(glPanel3d.window);
                int[] rowCol = adjustLimitRowCol(gather.superGatherRow, gather.superGatherCol);
                int currentRow = rowCol[0];
				int currentCol = rowCol[1];
                rowColIterator = getRowColIterator(1);
                while(rowColIterator.hasNext())
                {
                    rowCol = (int[])rowColIterator.next();
                    int row = rowCol[0];
                    int col = rowCol[1];
                    // hack until we figure out why row col is out of range
                    rowCol = adjustLimitRowCol(row, col);
                    row = rowCol[0];
                    col = rowCol[1];
                    getXYCoors(row, col, xyz);
                    if(row == currentRow && col == currentCol) continue;
                    if(velocityModel.hasVelocityProfile(row, col)) continue;
                    int numTracesInGather = getNTracesInGather(row, col);
                    if(numTracesInGather == 0)
                        continue;
                    else if (numTracesInGather < traceThreshold)
                    {
                        StsGLDraw.drawPoint(xyz, StsColor.BLACK, glPanel3d, size+buff, 1.0);
                        StsGLDraw.drawPoint(xyz, StsColor.RED, glPanel3d, size, 2.0);
                    }
                    else
                    {
                        StsGLDraw.drawPoint(xyz, StsColor.BLACK, glPanel3d, size+buff, 1.0);
                        StsGLDraw.drawPoint(xyz, StsColor.BLUE, glPanel3d, size, 2.0);
                    }
                }
                StsObjectRefList velocityProfiles = velocityModel.getVelocityProfiles();
                int nProfiles = velocityProfiles.getSize();
                for(int n = 0; n < nProfiles; n++)
                {
                    StsVelocityProfile profile = (StsVelocityProfile)velocityProfiles.getElement(n);
                    getXYCoors(profile.row, profile.col, xyz);
                    StsGLDraw.drawPoint(xyz, StsColor.BLACK, glPanel3d, size+buff, 1.0);
                    StsGLDraw.drawPoint(xyz, StsColor.GREEN, glPanel3d, size, 2.0);
                }
                getXYCoors(currentRow, currentCol, xyz);
                StsGLDraw.drawPoint(xyz, StsColor.BLACK, glPanel3d, size+buff, 2.0);
                StsGLDraw.drawPoint(xyz, StsColor.YELLOW, glPanel3d, size, 3.0);
            }
			catch (Exception e)
			{
				StsException.outputWarningException(this, "display", e);
			}
		}

		public void drawOnCursor2d(StsGLPanel3d glPanel3d, int dirNo, boolean axesFlipped)
		{
			if (dirNo != StsCursor3d.ZDIR) return;
            if(velocityModel == null) return;

			int size = lineSetClass.getDisplayPointSize();
			float[] xyz = new float[3];
			try
			{
				// Plot all the planned points
				StsSuperGather gather = getSuperGather(glPanel3d.window);
				int currentRow = gather.superGatherRow;
				int currentCol = gather.superGatherCol;
                rowColIterator = getRowColIterator(1);
                while(rowColIterator.hasNext())
                {
                    int[] rowCol = (int[])rowColIterator.next();
                    int row = rowCol[0];
                    int col = rowCol[1];
                    getXYCoors(row, col, xyz);
                    float x, y;
                    if (!axesFlipped)
                    {
                        x = xyz[0];
                        y = xyz[1];
                    }
                    else
                    {
                        x = xyz[1];
                        y = xyz[0];
                    }
                    if(row == currentRow && col == currentCol)
                    {
                        StsGLDraw.drawPoint2d(x, y, StsColor.BLACK, glPanel3d, size + 2, 0.0);
					    StsGLDraw.drawPoint2d(x, y, StsColor.YELLOW, glPanel3d, size, 2.0);
                        continue;
                    }
                    if(velocityModel.hasVelocityProfile(row, col))
                    {
                        StsGLDraw.drawPoint2d(x, y, StsColor.BLACK, glPanel3d, size + 2, 0.0);
						StsGLDraw.drawPoint2d(x, y, StsColor.GREEN, glPanel3d, size, 2.0);
                        continue;
                    }
                    int numTracesInGather = getNTracesInGather(row, col);
                    if (numTracesInGather < traceThreshold)
                    {
                        StsGLDraw.drawPoint2d(xyz[0], xyz[1], StsColor.BLACK, glPanel3d, size + 2, 0.0);
					    StsGLDraw.drawPoint2d(xyz[0], xyz[1], StsColor.RED, glPanel3d, size, 2.0);
                    }
                    else
                    {
                        StsGLDraw.drawPoint2d(xyz[0], xyz[1], StsColor.BLACK, glPanel3d, size + 2, 0.0);
                        StsGLDraw.drawPoint2d(xyz[0], xyz[1], StsColor.BLUE, glPanel3d, size, 2.0);
                    }
                }
			}
			catch (Exception e)
			{
				StsException.outputException("StsPreStackLineSet.VelStatus.display failed.", e, StsException.WARNING);
				return;
			}
		}

	}

    public void batchProcess(byte type)
    {
        Runnable runExport;
        Thread exportThread;

        switch(type)
        {
            case CVS:
                progressBarDialog = StsProgressBarDialog.constructor(currentModel.win3d, "Batch CVS Export", false);
                progressBarDialog.setLabelText("Exporting CVS Panels for PreStack3d: " + getName());
                progressBarDialog.pack();

                runExport = new Runnable()
                {
                    public void run()
                    {
                        batchProcessCVS();
                    }
                };
                exportThread = new Thread(runExport);
                exportThread.start();
                break;
            case SEMBLANCE:
                progressBarDialog = StsProgressBarDialog.constructor(currentModel.win3d, "Batch Semblance Export", false);
                progressBarDialog.setLabelText("Exporting semblance gathers for PreStack3d: " + getName());
                progressBarDialog.pack();

                runExport = new Runnable()
                {
                    public void run()
                    {
                        batchProcessSemblance();
                    }
                };
                exportThread = new Thread(runExport);
                exportThread.start();
                break;
            default:
                return;
        }
    }

    public StsBinaryFile openBatchOutputFile(byte type)
    {
        String dirname = currentModel.getProject().getDataFullDirString();
        String filename = "Unknown";
        try
        {
            switch(type)
            {
                case CVS:
                    filename = "cvs3d.bin." + this.getName();
                    if(this instanceof StsPreStackLineSet2d)
                        filename = "cvs2d.bin." + this.getName();
                    break;
                case SEMBLANCE:
                    filename = "semblance3d.bin." + this.getName();
                    if(this instanceof StsPreStackLineSet2d)
                        filename = "semblance2d.bin." + this.getName();
                    break;
                default:
                    return null;
            }
            StsBinaryFile bfile = new StsBinaryFile(StsFile.constructor(dirname, filename));
            bfile.openWrite();
            return bfile;
        }
        catch(Exception e)
        {
            StsException.outputException("StsPreStackLineSet.outputBatchHeaderFile failed.", e, StsException.WARNING);
            return null;
        }
    }

    public boolean outputBatchHeaderFile(byte type)
    {
        String filename = currentModel.getProject().getDataFullDirString();

        Properties props = new Properties(filterProperties, agcProperties, semblanceDisplayProperties,
                                          semblanceEditProperties, semblanceComputeProperties, semblanceRangeProperties,
                                          datumProperties, superGatherProperties, cvsProperties );
        try
        {
            switch(type)
            {
                case CVS:
                    filename = filename + "/cvs3d.txt." + this.getName();
                    if(this instanceof StsPreStackLineSet2d)
                        filename = filename + "/cvs2d.txt." + this.getName();
                    break;
                case SEMBLANCE:
                    filename = filename + "/semblance3d.txt." + this.getName();
                    if(this instanceof StsPreStackLineSet2d)
                        filename = filename + "/semblance2d.txt." + this.getName();
                    break;
                default:
                    return false;
            }
            StsParameterFile.writeObjectFields(filename, props);
        }
        catch(Exception e)
        {
            StsException.outputException("StsPreStackLineSet.outputBatchHeaderFile failed.", e, StsException.WARNING);
        }
        return true;
    }

    public void batchProcessCVS()
    {
        StsSuperGather gather = StsSuperGather.constructor(currentModel, this);
        outputBatchHeaderFile(CVS);
        StsBinaryFile outFile = openBatchOutputFile(CVS);
        try
        {
//            cvsProperties.initializeVelocity(this);
            progressBarDialog.setVisible(true);

            // Loop over every velocity analysis gather
            int row, col;
            int nLines = lines.length;
            for (int n = 0; n < nLines; n++)
            {
                StsPreStackLine line = lines[n];
                if (line.isInline)
                {
                    row = line.lineIndex;
                    int colMin = line.minGatherIndex;
                    int colMax = line.maxGatherIndex;
                    progressBarDialog.setProgressMax((colMax - colMin + 1));
                    int progress = 0;
                    for (col = colMin; col <= colMax; col++, progress++)
                    {
                        if(!gather.initializeSuperGather(row, col)) continue;
                        exportCVSGather(gather, outFile);
                        progressBarDialog.setProgress(progress);
                        progressBarDialog.setLabelText("Completed CVS for gather at CDP (" + (row+1) + "," + (col + 1) + ")");
                    }
                }
                else
                {
                    col = line.lineIndex;
                    int rowMin = line.minGatherIndex;
                    int rowMax = line.maxGatherIndex;
                    progressBarDialog.setProgressMax((rowMax - rowMin + 1));
                    int progress = 0;
                    for (row = rowMin; row <= rowMax; row++)
                    {
                        if(gather.initializeSuperGather(row, col))
                        {
                            exportCVSGather(gather, outFile);
                            progressBarDialog.setLabelText("Completed CVS for gather at CDP (" + (row+1) + "," + (col + 1) + ")");
                        }
                        progressBarDialog.setProgress(progress++);
                    }
                }
            }
            progressBarDialog.setLabelText("Completed CVS stacks.");
            progressBarDialog.dispose();
        }
        finally
        {
            // Close the output file
            outFile.close();
        }
    }

     public void exportCVSGather(StsSuperGather superGather, StsBinaryFile outFile)
    {
        float[][][] localCVSTraces = new float[cvsProperties.getNumberPanels()][][];

        boolean ignore = cvsProperties.getIgnoreSuperGather() || superGatherProperties.getIgnoreSuperGather();
        int tracesPerPanel = cvsProperties.getTracesPerPanel();
        if(!ignore)
            tracesPerPanel = superGather.getNumberGathers();

        superGather.checkInitializeSuperGatherGeometry();
        float velocity = semblanceRangeProperties.velocityMin;
        for (int k = 0; k < tracesPerPanel; k++)
        {
            if (!ignore)
            {
                localCVSTraces[k] = superGather.computeSuperGatherCVSTraces(nSlices, velocity);
            }
            else
            {
                localCVSTraces[k] = new float[tracesPerPanel][];
                localCVSTraces[k][0] = superGather.centerGather.computeConstantVelocityStackTrace(nSlices, velocity);
                for (int m = 1; m < tracesPerPanel; m++)
                {
                    localCVSTraces[k][m] = localCVSTraces[k][0];
                    outFile.setFloatValues(localCVSTraces[k][m], false);
                }
            }
            velocity = velocity + cvsProperties.getCvsVelocityStep(this);
        }
    }

    public void batchProcessSemblance()
    {
        StsSuperGather gather = StsSuperGather.constructor(currentModel, this);
        outputBatchHeaderFile(SEMBLANCE);
        StsBinaryFile outFile = openBatchOutputFile(SEMBLANCE);
        try
        {
//            cvsProperties.initializeVelocity(this);
            progressBarDialog.setVisible(true);

            // Loop over every velocity analysis gather
            int row, col;
            int nLines = lines.length;
            for (int n = 0; n < nLines; n++)
            {
                StsPreStackLine line = lines[n];
                if (line.isInline)
                {
                    row = line.lineIndex;
                    int colMin = line.minGatherIndex;
                    int colMax = line.maxGatherIndex;
                    progressBarDialog.setProgressMax((colMax - colMin + 1));
                    int progress = 0;
                    for (col = colMin; col <= colMax; col++)
                    {
                        if(gather.initializeSuperGather(row, col))
                        {
                            if(!exportSemblance(gather, outFile))
                                StsMessageFiles.infoMessage("Failed to compute semblance for CDP (" + (row+1) + "," + (col + 1) + ")");
                            else
                                progressBarDialog.setLabelText("Completed semblance for gather at CDP (" + (row+1) + "," + (col + 1) + ")");
                            progressBarDialog.setProgress(progress++);
                        }
                    }
                }
                else
                {
                    col = line.lineIndex;
                    int rowMin = line.minGatherIndex;
                    int rowMax = line.maxGatherIndex;
                    progressBarDialog.setProgressMax((rowMax - rowMin + 1));
                    int progress = 0;
                    for (row = rowMin; row <= rowMax; row++)
                    {
                        if(gather.initializeSuperGather(row, col))
                        {
                            if(!exportSemblance(gather, outFile))
                                StsMessageFiles.infoMessage("Failed to compute semblance for CDP (" + (row+1) + "," + (col + 1) + ")");
                            else
                                progressBarDialog.setLabelText("Completed Semblance for gather at CDP (" + (row+1) + "," + (col + 1) + ")");
                        }
                        progressBarDialog.setProgress(progress++);
                    }
                }
            }
            progressBarDialog.setLabelText("Completed Semblance Export.");
            progressBarDialog.dispose();
        }
        finally
        {
            outFile.close();
        }
    }

    public boolean exportSemblance(StsSuperGather superGather, StsBinaryFile outFile)
    {
        if (superGather.nSuperGatherTraces == 0)
        {
            superGather.createNullSemblanceBytes();
        }
        else
        {
            byte semblanceType = getSemblanceComputeProperties().semblanceType;
            if (semblanceType == StsSemblanceComputeProperties.SEMBLANCE_STANDARD)
            {
                if(!computeSemblanceStandard(superGather))
                    return false;
            }
            else if (semblanceType == StsSemblanceComputeProperties.SEMBLANCE_ENHANCED)
            {
                if(!computeSemblanceEnhanced(superGather))
                    return false;
            }
        }
        // Output the current semblance
        outFile.setByteValues(superGather.semblanceBytes);
        return true;
    }

    public boolean computeSemblanceStandard(StsSuperGather superGather)
    {
        return superGather.computeSemblanceStandardProcess(this, null);
    }

    public boolean computeSemblanceEnhanced(StsSuperGather gather)
    {
        return false;
    }

	 public boolean canExport()
	 {
		 return false;
		 // TODO When ready to finish off implementing the export here change to true.
		 //return true;
	 }

	 public boolean export()
	 {
		  try
		  {
				new StsProcessDismissDialog(currentModel.win3d, this, "Seismic Export", true);
				return true;
		  }
		  catch(Exception e)
		  {
				new StsMessage(currentModel.win3d, StsMessage.WARNING, "Seismic export failed.");
				return false;
		  }
	 }

	 public StsDialogFace getEditableCopy()
	 {
		  return (StsDialogFace)StsToolkit.copyObjectNonTransientFields(this);
	 }

	 public void dialogSelectionType(int type)
	 {
		  exportPanel.dialogSelectionType(type);
	 }

    public RowColIterator getRowColIterator(int direction)
    {
        return getRowColIterator(analysisRowStart, analysisColStart, direction);
    }

    public RowColIterator getRowColIterator(int row, int col, int direction)
    {
        if(rowColIterator == null)
            rowColIterator = new RowColIterator(row, col, direction);
        else
            rowColIterator.initialize(row, col, direction);
        return rowColIterator;
    }

    public class RowColIterator implements Iterator
    {
        public int row, col;
        int direction;
        int startRow, startCol;
        int endRow, endCol;
        boolean hasNext = false;

        RowColIterator(int row, int col, int direction)
        {
            initialize(row, col, direction);
        }

        public void initialize(int r, int c, int direction)
        {
            if(direction > 0)
            {
                row = StsMath.intervalRoundDown(r, analysisRowStart, analysisRowInc);
                col = StsMath.intervalRoundDown(c, analysisColStart, analysisColInc);
            }
            else
            {
                row = StsMath.intervalRoundUp(r, analysisRowStart, analysisRowInc);
                col = StsMath.intervalRoundUp(c, analysisColStart, analysisColInc);
            }
            endRow = nRows-1;
            endRow = StsMath.intervalRoundDown(endRow, analysisRowStart, analysisRowInc);
            endCol = getNColsForRow(row)-1;
            endCol = StsMath.intervalRoundDown(endCol, analysisColStart, analysisColInc);
            row = StsMath.minMax(this.row, analysisRowStart, endRow);
            col = StsMath.minMax(this.col, analysisColStart, endCol);
            startRow = row;
            startCol = col;
            this.direction = direction;
            hasNext = true;
        }

        public boolean hasNext()
        {
            return hasNext;
        }

        public Object next()
        {
            if(!hasNext) return null;
            int[] rowCol = new int[] { row, col };
            if(direction > 0)
            {
                col += analysisColInc;
                if(col > endCol)
                {
                    col = analysisColStart;
                    row += analysisRowInc;
                    if(row > endRow)
                        row = analysisRowStart;
                    endCol = getNColsForRow(row)-1;
                    endCol = StsMath.intervalRoundDown(endCol, analysisColStart, analysisColInc);
                }
            }
            else
            {
                col -= analysisColInc;
                if(col < analysisColStart)
                {
                    row -= analysisRowInc;
                    if(row < analysisRowStart)
                        row = endRow;
                    endCol = getNColsForRow(row)-1;
                    endCol = StsMath.intervalRoundDown(endCol, analysisColStart, analysisColInc);
                    col = endCol;
                }
            }
            hasNext = row != startRow || col != startCol;
            return rowCol;
        }

        public void remove()
        {
        }
    }

    class Properties
    {
        StsFilterProperties filterProperties = null;
        StsAGCPreStackProperties agcProperties = null;
        StsSemblanceDisplayProperties semblanceDisplayProperties = null;
        StsSemblanceDisplayProperties semblanceEditProperties = null;
        StsSemblanceComputeProperties semblanceComputeProperties = null;
        StsSemblanceRangeProperties semblanceRangeProperties = null;
        StsDatumProperties datumProperties = null;
        StsSuperGatherProperties superGatherProperties = null;
        StsCVSProperties cvsProperties = null;

        Properties(StsFilterProperties filterProperties, StsAGCPreStackProperties agcProperties,
                   StsSemblanceDisplayProperties semblanceDisplayProperties, StsSemblanceDisplayProperties semblanceEditProperties,
                   StsSemblanceComputeProperties semblanceComputeProperties, StsSemblanceRangeProperties semblanceRangeProperties,
                   StsDatumProperties datumProperties, StsSuperGatherProperties superGatherProperties,
                   StsCVSProperties cvsProperties )
        {
            this.filterProperties = filterProperties;
            this.agcProperties = agcProperties;
            this.semblanceDisplayProperties = semblanceDisplayProperties;
            this.semblanceEditProperties = semblanceEditProperties;
            this.semblanceComputeProperties = semblanceComputeProperties;
            this.semblanceRangeProperties = semblanceRangeProperties;
            this.datumProperties = datumProperties;
            this.superGatherProperties = superGatherProperties;
            this.cvsProperties = cvsProperties;
        }
	}
    /** Initialization before data is fetched from source interface. Aadded for compatibility with StsDataLineSetFace interface */
    public void initializeDataSet()
    {
    }
    
    /**
     * get CDP number row and column number
     * @param row sequential counter of row
     * @param col sequential counter of the gather along row
     * @return CDP number or -1 if line can't be found at the row and column number
     */
    public int getCDP(int row, int col)
    {
        StsPreStackLine line = this.getDataLine(row, col);
        if (line == null) return -1;
        return line.getCDP(col);
    }
    public void setSuperGatherToCursor(StsWin3dBase window, StsCursor3d cursor3d)
    {
        float[] xyz = cursor3d.getCurrentCoordinates();
        int[] rowCol = getRowColFromCoors(xyz[0], xyz[1]);
        setCurrentDataRowCol(window, rowCol[0], rowCol[1]);
    }
    public int getCDPFromRowNumColNum(int rowNum, int colNum)
    {
        int row = this.getRowFromRowNum(rowNum);
        int col = this.getColFromColNum(colNum);
        StsPreStackLine line = this.getDataLine(row, col);
        if (line == null) return -1;
        double[] xlines = null;
        double[] cdps = null;
        try
        {
            xlines = line.getAttributeArray(StsSEGYFormat.XLINE_NO);
            cdps = line.getAttributeArray(StsSEGYFormat.CDP);
        } catch (FileNotFoundException e)
        {
            e.printStackTrace();
        }
        if (xlines == null) return -1;
        for (int i=0; i<xlines.length; i++) 
        {
            if (colNum == xlines[i]) return (int)cdps[i];
        }
        return -1;
    }
    public StsWiggleDisplayProperties getPanelWiggleDisplayProperties()
    {
        return lineSetClass.getPanelWiggleDisplayProperties();
    }
}
