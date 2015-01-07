package com.Sts.DBTypes;

/**
 * <p>Title: Workflow development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2001</p>
 * <p>Company: 4D Systems LLC</p>
 * @author unascribed
 * @version 1.0
 */

import com.Sts.Actions.Wizards.Color.*;
import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.text.*;
import java.util.*;

public class StsVsp extends StsSeismicBoundingBox implements ItemListener, StsTreeObjectI, StsTextureSurfaceFace //, StsEfficientRepaintable
{
	// these members are persistent, but not loaded from seis3d.txt.name file
	protected StsColorscale colorscale;
	protected StsWell well = null;

	transient boolean readoutEnabled = false;
	transient StsFile out = null;
	transient OutputStream os = null;
	transient BufferedOutputStream bos = null;
	transient ByteArrayOutputStream baos = null;
	transient DataOutputStream ds = null;
//	transient byte exportType = StsSeismicExportPanel.BYTE;
//	transient byte exportNullType = StsSeismicExportPanel.NULL_NONE;
//	transient float exportScale = 1.0f;
//	transient String exportTextHeader = null;
	//transient public StsVspView vspView;

	transient StsProgressBarDialog progressBarDialog = null;

	// the following are initialized after reading parameters file or reloading database
	transient protected StsSpectrumDialog spectrumDialog = null;
	transient protected boolean spectrumDisplayed = true;

	transient protected int colorListNum = 0;
	transient protected boolean colorListChanged = true;

	transient int sampleSize = 1;
    transient boolean attributeInMdepth = false;
    transient byte attributeDomain = StsParameters.TD_DEPTH;

	transient Vector itemListeners = null;
	transient StsActionListeners actionListeners = null;

	transient StsTimer timer = null;
	transient boolean runTimer = false;

	transient static final boolean debug = false;
	transient static final boolean voxelDebug = false;
	transient static final boolean wiggleDebug = false;

	/** Total samples in each of 255 steps */
	transient private int dataCnt[] = new int[255];
	transient private int ttlHistogramSamples = 0;

	transient int displayListNum = 0;

	transient float[] floatData = null; // data read for this vsp
	transient float[] scaledFloatData = null; // scaled, agc'd etc trace data
	transient byte[] byteData = null; // byte data converted from fdata

	private transient StsTextureTiles textureTiles = null;
	transient private boolean textureChanged = false;
	transient float azimuth;
	transient int panelWidth = 100;
	transient float[][] panelRanges;
	transient boolean isDisplayWiggles = false;
	transient boolean isPixelMode = false;
	transient String exportDirectory = null;
	transient StsVspClass vspClass;
	transient public StsColorList seismicColorList = null;

	public final static byte nullByte = StsParameters.nullByte;
	public final static int nullUnsignedInt = 255;

	static protected StsObjectPanel objectPanel = null;

	protected boolean displayAxis = true;
	protected boolean displayVSPs = true;
	protected boolean displayVAR = true;
	protected boolean displayWiggles = false;
	protected boolean contourColors = true;
	protected int wiggleToPixelRatio = 1;
	//protected float tracesPerInch = 25.f;
	protected transient float currentTracesPerInch = 25.f;
    protected transient float currentInchesPerSecond = 3.f;

//	private StsWiggleDisplayProperties wiggleProperties = null;
	public StsFilterProperties filterProperties = null;
	public StsAGCProperties agcProperties = null;
    public StsDatumProperties datumProperties = null;

	transient public String attributeFilename = null;

	final static public String attributeFilePrefix = "seis.bin.";

//	static public final StsFloatFieldBean tpiBean = new StsFloatFieldBean(StsVsp.class, "tracesPerInch", 1, 250, "Traces per Inch:");

    static public StsFieldBean[] displayFields = null;
    static StsDateFieldBean bornField = new StsDateFieldBean(StsVsp.class, "bornDate", "Born Date:");
    static StsDateFieldBean deathField = new StsDateFieldBean(StsVsp.class, "deathDate", "Death Date:");

    static StsEditableColorscaleFieldBean colorscaleBean;
    static public StsFieldBean[] propertyFields = null;

    public String getGroupname()
    {
        return groupVsp;
    }

    public StsVsp()
	{
	}

	public StsVsp(boolean persistent)
	{
		super(persistent);
	}

	public StsVsp(String[] args)
	{
		setYInc(Float.parseFloat(args[0]));
		setXInc(Float.parseFloat(args[1]));
		setAngle(Float.parseFloat(args[2]));
//        volume.setIsXLineCCW(true);
		initialize((StsModel)null);
	}

	private StsVsp(String directory, String filename, StsModel model) throws FileNotFoundException, StsException
	{
		super(false);
		String pathname = directory + filename;
		if(!(new File(pathname)).exists())
		{
			throw new FileNotFoundException();
		}
		StsParameterFile.initialReadObjectFields(pathname, this, StsSeismicBoundingBox.class, StsMainObject.class);
		byte zDomainByte = StsParameters.getZDomainFromString(zDomain);
		/* SAJ
		   if(zDomainByte != StsParameters.TD_TIME)
		   {
		 new StsMessage(currentModel.win3d, StsMessage.WARNING, "VSP file " + filename + " is not in time. Not supported.");
		 throw new StsException(StsException.WARNING, "VSP volume not in time.");
		   }
		 */
		//model.getProject().addToProject(this);
		setName(getStemname());
		initializeColorscale();
		stsDirectory = directory;
		vspClass = (StsVspClass)currentModel.getStsClass(StsVsp.class);
//		setWiggleProperties(new StsWiggleDisplayProperties(this, vspClass.defaultWiggleProperties, "wiggleProperties"));
		if(!initialize(model))
		{
			throw new FileNotFoundException(pathname);
		}
		isVisible = true;
		this.checkDataFile();

		filterProperties = new StsFilterProperties(this, vspClass.defaultFilterProperties, "filterProperties");
		agcProperties = new StsAGCProperties(this, vspClass.defaultAGCProperties, "agcProperties");
        datumProperties = new StsDatumProperties(this, vspClass.defaultDatumProperties, "datumProperties");
		addToModel();
		//currentModel.add(wiggleProperties);
	}

	static public StsVsp constructor(StsFile file, StsModel model)
	{
		return StsVsp.constructor(file.getDirectory(), file.getFilename(), model);
	}

	static public StsVsp constructor(String directory, String filename, StsModel model)
	{
		StsVsp volume = null;
		try
		{
			volume = new StsVsp(directory, filename, model);
			model.setCurrentObject(volume);
//			volume.refreshObjectPanel();
			return volume;
		}
		catch(FileNotFoundException e)
		{
			return null;
		}
		catch(Exception e)
		{
			StsException.outputException("StsVsp.constructor() failed.", e, StsException.WARNING);
			return null;
		}
	}

	public boolean initialize(StsFile file, StsModel model)
	{
		try
		{
			String pathname = file.getDirectory() + file.getFilename();
			StsParameterFile.initialReadObjectFields(pathname, this, StsSeismicBoundingBox.class, StsBoundingBox.class);
			setName(getStemname());
			initializeColorscale();
			stsDirectory = file.getDirectory();
			if(!initialize(model))
			{
				return false;
			}
			isVisible = true;

			// take *current* class settings as defaults
			{
				setDisplayVSPs(vspClass.getDisplayVSPs());
				setDisplayAxis(vspClass.getDisplayAxis());
				setIsPixelMode(vspClass.getIsPixelMode());
				setDisplayVAR(vspClass.getDisplayVAR());
				setContourColors(vspClass.getContourColors());
				//setDisplayWiggles(vspClass.getDisplayWiggles());
				setWiggleToPixelRatio(vspClass.getWiggleToPixelRatio());

			/* following not necessary since all wiggleProperties are initialized to current vspClass wiggleProperties in constructor.
			   TJL 8/30/06.
				setTracesPerInch(vspClass.getTracesPerInch());
                setInchesPerSecond(vspClass.getInchesPerSecond());
	        */
			}
			return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsVsp.loadFile() failed.",
										 e, StsException.WARNING);
			return false;
		}
	}

    // if crossline+ direction is 90 degrees CCW from inline+, this is isColCCW; otherwise not
	// angle is from X+ direction to inline+ direction (0 to 360 degrees)
	public boolean initialize(StsModel model)
	{
		try
		{
			super.initialize(currentModel);
			initializeColorscale();
            initializePropertyPanels();
			setDataHistogram();
			vspClass = (StsVspClass)currentModel.getStsClass(StsVsp.class);
			wiggleToPixelRatio = vspClass.getWiggleToPixelRatio();
			byte zDomainByte = StsParameters.getZDomainFromString(zDomain);
			zDomain = StsParameters.TD_ALL_STRINGS[zDomainByte];
            currentTracesPerInch = getTracesPerInch();
            currentInchesPerSecond = getInchesPerSecond();
			return currentModel.getProject().checkSetZDomain(zDomainByte, zDomainByte);

		}
		catch(Exception e)
		{
			StsException.outputException("StsVsp.classInitialize(model) failed.", e, StsException.WARNING);
//			StsMessage.printMessage("Failed to find file. Error: " + e.getMessage());
			return false;
		}
	}

	private boolean initialize()
	{
		super.initialize(currentModel);
		setDataHistogram();
		byte zDomainByte = StsParameters.getZDomainFromString(zDomain);
		zDomain = StsParameters.TD_ALL_STRINGS[zDomainByte];
		vspClass = (StsVspClass)currentModel.getStsClass(StsVsp.class);
		return currentModel.getProject().checkSetZDomain(zDomainByte, zDomainByte);
	}

    private void initializePropertyPanels()
    {
        if(filterProperties != null) filterProperties.setParentObject(this);
        if(agcProperties != null) agcProperties.setParentObject(this);
        if(datumProperties != null) datumProperties.setParentObject(this);
        if(datumProperties != null) datumProperties.buildBeanLists(this);
    }

    public void displayWiggleProperties()
	{
		getWiggleDisplayProperties().displayWiggleProperties("Edit " + getName() + " Wiggle Properties");
	}

	public void displayAGCProperties()
	{
		displayWiggleProperties();
	}

	/**
	 * Set the pixel mode for VSP displays. When set to off, pixels will be interpolated between samples.
	 * When set to on, there will be hard edges between samples.
	 */

	public void setIsPixelMode(Boolean b)
	{
		setIsPixelMode(b.booleanValue());
	}

	public void setIsPixelMode(boolean isPixelMode)
	{
		if(this.isPixelMode == isPixelMode)return;
		this.isPixelMode = isPixelMode;
		dbFieldChanged("isPixelMode", isPixelMode);
		currentModel.viewObjectChangedAndRepaint(this, this);

	}

	public boolean getIsPixelMode()
	{
		return isPixelMode;
	}

	/**
	 * Set the display axis toggle
	 */
	public void setDisplayAxis(Boolean b)
	{
		setDisplayAxis(b.booleanValue());
	}

	public void setDisplayAxis(boolean displayAxis)
	{
		if(this.displayAxis == displayAxis)return;
		this.displayAxis = displayAxis;
		dbFieldChanged("displayAxis", displayAxis);

		currentModel.viewObjectChangedAndRepaint(this, this);
	}

	public boolean getDisplayAxis()
	{
		return displayAxis;
	}

	public boolean getDisplayVSPs()
	{
		return displayVSPs;
	}

	public void setDisplayVSPs(Boolean b)
	{
		setDisplayVSPs(b.booleanValue());
	}

	public void setDisplayVSPs(boolean displayVSPs)
	{
		if(this.displayVSPs == displayVSPs)return;
		this.displayVSPs = displayVSPs;
		dbFieldChanged("displayVSPs", displayVSPs);

		currentModel.viewObjectChangedAndRepaint(this, this);
	}

    public boolean getUseShader() { return vspClass.getContourColors(); }
    public int getDefaultShader() { return StsJOGLShader.ARB_TLUT_NO_SPECULAR_LIGHTS; }
	public void setContourColors(Boolean b)
	{
		setContourColors(b.booleanValue());
	}

	public void setContourColors(boolean contourColors)
	{
		if(this.contourColors == contourColors)return;
		this.contourColors = contourColors;
		dbFieldChanged("contourColors", contourColors);

		currentModel.viewObjectChangedAndRepaint(this, this);
	}

	public boolean getContourColors()
	{
		return contourColors;
	}
    /** display variable area toggle */
	public void setDisplayVAR(Boolean b)
	{
		setDisplayVAR(b.booleanValue());
	}

	public void setDisplayVAR(boolean displayVAR)
	{
		if(this.displayVAR == displayVAR)return;
		this.displayVAR = displayVAR;
		dbFieldChanged("displayVAR", displayVAR);

		currentModel.viewObjectChangedAndRepaint(this, this);
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
	public void setWiggleToPixelRatio(Integer wiggleToPixelRatio)
	{
		setWiggleToPixelRatio(wiggleToPixelRatio.intValue());
	}

	public void setWiggleToPixelRatio(int wiggleToPixelRatio)
	{
		if(this.wiggleToPixelRatio == wiggleToPixelRatio)
			return;
		this.wiggleToPixelRatio = wiggleToPixelRatio;
		dbFieldChanged("wiggleToPixelRatio", wiggleToPixelRatio);
		//setDisplayField("wiggleToPixelRatio", wiggleToPixelRatio);

		currentModel.viewObjectChangedAndRepaint(this, this);
	}

	public int getWiggleToPixelRatio()
	{
		return wiggleToPixelRatio;
	}
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
		if(!setval)
		{
			tpiBean.setValueObjectFromString(new Float(tracesPerInch).toString());

			return;
		}
		if(this.tracesPerInch == tracesPerInch)return;
		this.tracesPerInch = tracesPerInch;
		dbFieldChanged("tracesPerInch", tracesPerInch);

		tpiBean.setValueObjectFromString(new Float(tracesPerInch).toString());

		currentModel.viewObjectChanged(this);
		currentModel.viewObjectRepaint(this);

	}

	public float getTracesPerInch()
	{
		return tracesPerInch;
	}

	 */

    public void setInchesPerSecond(Float inchesPerSecond)
    {
        getWiggleDisplayProperties().setInchesPerSecond(inchesPerSecond.floatValue());
    }

    public void setInchesPerSecond(float inchesPerSecond)
    {
        getWiggleDisplayProperties().setInchesPerSecond(inchesPerSecond);
    }

    public float getInchesPerSecond()
    {
        return getWiggleDisplayProperties().getInchesPerSecond();
    }

    public float getCurrentInchesPerSecond()
    {
        return currentInchesPerSecond;
    }

    public void setTracesPerInch(Float tracesPerInch)
    {
        getWiggleDisplayProperties().setTracesPerInch(tracesPerInch.floatValue());
    }

    public void setTracesPerInch(float tracesPerInch)
    {
        getWiggleDisplayProperties().setTracesPerInch(tracesPerInch);
    }

    public void setTracesPerInch(float tracesPerInch, boolean setval)
    {
        getWiggleDisplayProperties().setTracesPerInch( tracesPerInch);

        currentModel.viewObjectChangedAndRepaint(this, this);
    }

    public float getTracesPerInch()
    {
        return getWiggleDisplayProperties().getTracesPerInch();
	}

	public float getCurrentTracesPerInch()
	{
		return currentTracesPerInch;
	}

	public String getSegyFileDate()
	{
		if(segyLastModified == 0)
		{
			File segyFile = new File(segyDirectory + segyFilename);
			if(segyFile != null)
			{
				segyLastModified = segyFile.lastModified();
			}
		}
		DateFormat dateFormat = DateFormat.getDateTimeInstance();
		return dateFormat.format(new Date(segyLastModified));
	}

	public String getDate()
	{
		return null;
	}

	public String getLabel()
	{
		return stemname;
	}

	public int getIntValue(int row, int col, int slice)
	{
		if(row < 0 || row >= nRows || col < 0 || col >= nCols || slice < 0 ||
		   slice >= nSlices)
		{
			StsException.systemError("StsVsp.getValue() failed for row: " + row + " col: " + col + " plane: " + slice);
			return 0;
		}

		byte[] rowData = (byte[])readPlaneData(StsCursor3d.YDIR, getYCoor(row));
		if(rowData == null)
		{
			return 0;
		}
		byte signedByteValue = rowData[col * nSlices + slice];
		return StsMath.signedByteToUnsignedInt(signedByteValue);
	}

	/*
	 public float getValue(float x, float y, float z)
	 {
	  int row = getNearestBoundedRowCoor(y);
	  int col = getNearestBoundedColCoor(x);
	  int plane = getNearestBoundedSliceCoor(z);
	  float f = (float)intValue(row, col, plane);
	  return dataMin + (f / 254) * (dataMax - dataMin);
	 }

	 public float getValue(int row, int col, int slice)
	 {
	  float f = (float)intValue(row, col, slice);
	  return dataMin + (f / 254) * (dataMax - dataMin);
	 }

	 public int intValue(float x, float y, float z)
	 {
	  int row = getNearestBoundedRowCoor(y);
	  int col = getNearestBoundedColCoor(x);
	  int plane = getNearestBoundedSliceCoor(z);
	  return intValue(row, col, plane);
	 }

	 public int intValue(float[] xyz)
	 {
	  int row = getNearestBoundedRowCoor(xyz[1]);
	  int col = getNearestBoundedColCoor(xyz[0]);
	  int plane = getNearestBoundedSliceCoor(xyz[2]);
	  return intValue(row, col, plane);
	 }

	 public int intValue(int row, int col, float z)
	 {
	  int v0, v1, v;

	  float kF = (float)((z - zMin) / zInc);
	  int k = (int)kF;
	  float dk = kF - k;
	  if(dk < 0.5f)
	  {
	   return intValue(row, col, k);
	  }
	  else
	  {
	   return intValue(row, col, k + 1);
	  }
	 }
	 */
	static public void close(StsVsp volume)
	{
		if(volume == null)return;
		volume.close();
	}

	public void close()
	{
	}

	public String getSegyFilename()
	{
		return segyFilename;
	}

	public void setSegyFilename(String segyFilename)
	{
		this.segyFilename = segyFilename;
	}

	public String getSegyDirectory()
	{
		return segyDirectory;
	}

	public void setSegyDirectory(String segyDirectory)
	{
		this.segyDirectory = segyDirectory;
	}

	public String getStsDirectory()
	{
		return stsDirectory;
	}

	public void setStsDirectory(String stsDirectory)
	{
		this.stsDirectory = stsDirectory;
	}

	public void setStemname(String stemname)
	{
		this.stemname = stemname;
	}

	public boolean getIsXLineCCW()
	{
		return isXLineCCW;
	}

	public void setIsXLineCCW(boolean isXLineCCW)
	{
		this.isXLineCCW = isXLineCCW;
	}

	public boolean getIsRegular()
	{
		return isRegular;
	}

	public void setIsRegular(boolean isRegular)
	{
		this.isRegular = isRegular;
	}

//    public float getSampleSpacing() { return getZInc(); }
//    public void setSampleSpacing(float sampleSpacing) { setZInc(sampleSpacing); }

//    public void setLineAngle(float lineAngle) { setAngle(lineAngle); }
//    public float getLineAngle() { return getAngle(); }

	public StsSpectrum getSpectrum()
	{
		return colorscale.getSpectrum();
	}

	public StsColorscale getColorscale()
	{
		setDataHistogram();
		return colorscale;
	}

	public StsColorscale getColorscaleWithName(String name)
	{
		if(name.equals(colorscale.getName()))
			return colorscale;
		else
			return null;
	}

	public void setDataHistogram()
	{
		if(colorscaleBean == null) return;
		colorscaleBean.setHistogram(dataHist);
	}

	/*
	 public void setColorscale(StsColorscale colorscale)
	 {
	  fieldChanged("colorscale", colorscale);
	  currentModel.win3dDisplayAll();
	 }
	 */
	public void itemStateChanged(ItemEvent e)
	{
		if(e.getItem() instanceof StsColorscale)
		{
			int id = e.getID();
			currentModel.displayIfCursor3dObjectChanged(this);
			fireItemStateChanged(e);
		}
		else
		{
			ItemEvent event = null;
			fireItemStateChanged(event);
		}
		return;
	}

	public void actionPerformed(ActionEvent e)
	{
		if(e.getSource() instanceof StsColorscale)
		{
			seismicColorList.setColorListChanged(true);
			currentModel.win3dDisplayAll();
			fireActionPerformed(e);
		}
		else
			fireItemStateChanged(null);
		return;
	}

	public synchronized void addActionListener(ActionListener listener)
	{
		if(actionListeners == null)actionListeners = new StsActionListeners();
		actionListeners.add(listener);
	}

	public synchronized void removeActionListener(ActionListener listener)
	{
		if(actionListeners == null)return;
		actionListeners.remove(listener);
	}

	protected void fireActionPerformed(ActionEvent e)
	{
		if(actionListeners == null)return;
		actionListeners.fireActionPerformed(e);
	}

	public void initializeColorscale()
	{
		try
		{
			if(colorscale == null)
			{
				StsSpectrumClass spectrumClass = currentModel.getSpectrumClass();
				colorscale = new StsColorscale("VSP", spectrumClass.getSpectrum(StsSpectrumClass.SPECTRUM_RWB), dataMin, dataMax);
				colorscale.setEditRange(dataMin, dataMax);
			}
			seismicColorList = new StsColorList(colorscale);
			colorscale.addActionListener(this);
		}
		catch(Exception e)
		{
			StsException.outputException("StsSeismicVolume.initializeColorscale() failed.", e, StsException.WARNING);
		}
	}

	public void setSpectrumDialog(StsSpectrumDialog spectrumDialog)
	{
		this.spectrumDialog = spectrumDialog;
	}

	public boolean getIsVisibleOnCursor()
	{
		return isVisible;
	}

	public boolean getReadoutEnabled()
	{
		return readoutEnabled;
	}

	public void setReadoutEnabled(boolean enabled)
	{
		readoutEnabled = enabled;
		return;
	}

	/*
	 public final Color getColor(float[] xyz)
	 {
	  int index = intValue(xyz);
	  return getColor(index);
	 }
	 */
	public final Color getColor(int colorIndex)
	{
		return colorscale.getColor(colorIndex);
	}

	/*
	 public final Color getColor(int row, int col, float z)
	 {
	  int index = intValue(row, col, z);
	  return getColor(index);
	 }

	 public final Color getColor(float x, float y, float z)
	 {
	  int index = intValue(x, y, z);
	  return getColor(index);
	 }

	 public StsColor getStsColor(float rowF, float colF, int plane)
	 {
	  int row = StsMath.roundOffInteger(rowF);
	  int col = StsMath.roundOffInteger(colF);
	  int v = intValue(row, col, plane);
	  return colorscale.getStsColor(intValue(row, col, plane));
	 }
	 */
	public void seismicChanged()
	{
		colorListChanged = true;
	}

	public static void main(String[] args)
	{
		try
		{
			String pathname = System.getProperty("user.dirNo");
            StsFileChooser chooser = StsFileChooser.createFileChooserPrefix(null, "select a seismic volume", pathname, "seis3d.txt");
			chooser.show();
			File[] files = chooser.getFiles();
			if(files == null || files[0] == null)
			{
				System.out.println("No file selected.");
				System.exit(0);
			}
			File file = files[0];
			pathname = file.getPath();
			StsVsp volume = new StsVsp();
			StsParameterFile.initialReadObjectFields(pathname, volume, StsSeismicBoundingBox.class, StsBoundingBox.class);
		}
		catch(Exception e)
		{
			StsException.outputException("StsVsp.main() failed.", e,
										 StsException.WARNING);
		}
		finally
		{
			System.exit(0);
		}
	}

	public void setIsVisible(boolean isVisible)
	{
		this.isVisible = isVisible;

		currentModel.viewObjectChangedAndRepaint(this, this);
	}

	public boolean getIsVisible()
	{
		return isVisible;
	}

//    public StsRotatedGridBoundingBox getBoundingBox() { return new StsRotatedGridBoundingBox(this); }

	public byte[] getCurrentCursorPlaneData()
	{
		return null;
	}

	public StsFieldBean[] getDisplayFields()
	{
        if(displayFields != null)
        {
			bornField = new StsDateFieldBean(StsVsp.class, "bornDate", "Born Date:");
			deathField = new StsDateFieldBean(StsVsp.class, "deathDate", "Death Date:");
            displayFields = new StsFieldBean[]
		    {
                new StsBooleanFieldBean(StsVsp.class, "displayVSPs", "Display on Well"),
				bornField, deathField,
                new StsBooleanFieldBean(StsVsp.class, "displayAxis", "Plot Axis"),
                new StsBooleanFieldBean(StsVsp.class, "isPixelMode", "Pixel Display Mode"),
                new StsBooleanFieldBean(StsVsp.class, "displayVAR", "Variable Area"),
                new StsBooleanFieldBean(StsVsp.class, "contourColors", "Contoured Seismic Colors"),
                // moved to wiggle props new StsBooleanFieldBean(StsVsp.class, "displayWiggles", "Wiggle Traces"),
                new StsIntFieldBean(StsVsp.class, "wiggleToPixelRatio", 1, 100, "Wiggle to Pixel Ratio:"),
                //tpiBean
	        };
        }
        return displayFields;
	}

	public StsFieldBean[] getPropertyFields()
	{
        if(propertyFields == null)
        {
	        colorscaleBean = new StsEditableColorscaleFieldBean(StsVsp.class, "colorscale");
            propertyFields = new StsFieldBean[]
            {
                new StsStringFieldBean(StsVsp.class, "name", true, "Name"),
        //        new StsStringFieldBean(StsVsp.class, "owner", false, "Owner:"),
        //        new StsStringFieldBean(StsVsp.class, "field", false, "Field:"),
        //        new StsStringFieldBean(StsVsp.class, "acquisitionDate", false, "Acquisition Date:"),
        //        new StsStringFieldBean(StsVsp.class, "acquiredBy", false, "Acquired By:"),
        //        new StsStringFieldBean(StsVsp.class, "processor", false, "Processor:"),
        //        new StsStringFieldBean(StsVsp.class, "processDate", false, "Process Date:"),

                new StsStringFieldBean(StsVsp.class, "zDomainString", false, "Z Domain"),
                new StsStringFieldBean(StsVsp.class, "segyFilename", false, "SEGY Filename"),
                new StsStringFieldBean(StsVsp.class, "segyFileDate", false, "SEGY creation date"),
                new StsIntFieldBean(StsVsp.class, "nCols", false, "Number of Crosslines"),
                new StsIntFieldBean(StsVsp.class, "nCroppedSlices", false, "Number of Samples"),
                new StsDoubleFieldBean(StsVsp.class, "xOrigin", false, "X Origin"),
                new StsDoubleFieldBean(StsVsp.class, "yOrigin", false, "Y Origin"),
                new StsFloatFieldBean(StsVsp.class, "xInc", false, "X Inc"),
                new StsFloatFieldBean(StsVsp.class, "zInc", false, "Z Inc"),
                new StsFloatFieldBean(StsVsp.class, "dataMin", false, "Data Min"),
                new StsFloatFieldBean(StsVsp.class, "dataMax", false, "Data Max"),
                new StsFloatFieldBean(StsVsp.class, "dataAvg", false, "Data Avg"),
                colorscaleBean,
                new StsButtonFieldBean("VSP Properties", "Edit VSP properties.", StsVsp.class, "displayVspPropertiesDialog")
            };
        }
        return propertyFields;
	}

	public Object[] getChildren()
	{
		return new Object[0];
	}

	public StsObjectPanel getObjectPanel()
	{
		if(objectPanel == null)
		{
			objectPanel = StsObjectPanel.constructor(this, true);
		}
		return objectPanel;
	}

	public void treeObjectSelected()
	{
		currentModel.getStsClass(StsVsp.class).selected(this);
		currentModel.getGlPanel3d().checkAddView(StsView3d.class);

		currentModel.viewObjectChangedAndRepaint(this, this);
	}

	public boolean anyDependencies()
	{
		StsCrossplot[] cp = (StsCrossplot[])currentModel.getCastObjectList(
			StsCrossplot.class);
		for(int n = 0; n < cp.length; n++)
		{
			StsSeismicBoundingBox[] volumes = cp[n].getVolumes();
			for(int j = 0; j < cp[n].volumes.getSize(); j++)
			{
				if(this == cp[n].volumes.getElement(j))
				{
					StsMessageFiles.infoMessage("Seismic PostStack3d " + getName() +
												" used by Crossplot " +
												cp[n].getName());
					return true;
				}
			}
		}
		return false;
	}

	public boolean delete()
	{
		close();
		// the class could be either StsVsp or StsVirtualVolume
		StsClass stsClass = currentModel.getStsClass(getClass());
		StsObject currentObject = stsClass.getCurrentObject();
		stsClass.delete(this);
		super.delete(); // sets this.index to -1

		if(currentObject == this)
		{
			int nVolumes = stsClass.getSize();
			if(nVolumes > 0)
			{
				currentObject = stsClass.getElement(nVolumes - 1);
			}
			else
			{
				currentObject = null;

			}
			stsClass.setCurrentObject(currentObject);
			if(currentModel.getGlPanel3d() != null)
				currentModel.toggleOnCursor3dObject(currentObject);
		}
		return true;
	}

	public void setColorscale(StsColorscale colorscale)
	{
		seismicColorList.setColorscale(colorscale);
		this.colorscale = colorscale;

		currentModel.viewObjectChangedAndRepaint(this, this);
	}

	public boolean setGLColorList(GL gl, boolean nullsFilled, int shader)
	{
		return seismicColorList.setGLColorList(gl, nullsFilled, shader);
	}

	/** Gets the display type (Cell or Grid); */
//	public byte getDisplayType() { return StsSurfaceDisplayable.ORTHO_GRID_CENTERED; }
	/** This property can be isVisible on this surace. */
	public boolean isDisplayable()
	{
		return true;
	}

	/** Call for getting a cell-centered color */
	public Color getCellColor(int row, int col, int layer)
	{
		return Color.RED;
	}

	/** Call for getting a grid-centered color */
	/*
	 public Color getGridColor(int row, int col, float z)
	 {
	  return getColor(row, col, z);
	 }
	 */
	public Color[][] get2dColorArray(StsRotatedGridBoundingBox surfaceBoundingBox, float[][] z, float zOffset)
	{
		int nSurfaceRows = surfaceBoundingBox.nRows;
		int nSurfaceCols = surfaceBoundingBox.nCols;

		Color[][] colors = new Color[nSurfaceRows][nSurfaceCols];
		Color[] colorscaleColors = colorscale.getNewColorsInclTransparency();
		int value;
		int currentK = -1;
		byte[] planeData = null;

		if(timer == null)
		{
			timer = new StsTimer();
		}
		timer.start();

		int surfaceRowMin = surfaceBoundingBox.getNearestBoundedRowCoor(yMin);
		int surfaceRowMax = surfaceBoundingBox.getNearestBoundedRowCoor(yMax);
		int surfaceColMin = surfaceBoundingBox.getNearestBoundedColCoor(xMin);
		int surfaceColMax = surfaceBoundingBox.getNearestBoundedColCoor(xMax);

		Color nullColor = colorscaleColors[0];

		for(int surfaceRow = 0; surfaceRow < nSurfaceRows; surfaceRow++)
		{
			if(surfaceRow < surfaceRowMin || surfaceRow > surfaceRowMax)
			{
				for(int surfaceCol = 0; surfaceCol < nSurfaceCols; surfaceCol++)
				{
					colors[surfaceRow][surfaceCol] = nullColor;
				}
			}
			else
			{
				int row = surfaceRow - surfaceRowMin;
				for(int surfaceCol = 0; surfaceCol < nSurfaceCols; surfaceCol++)
				{
					if(surfaceCol < surfaceColMin || surfaceCol > surfaceColMax)
					{
						colors[surfaceRow][surfaceCol] = nullColor;
					}
					else
					{
						int col = surfaceCol - surfaceColMin;

						float kF = (float)((z[surfaceRow][surfaceCol] + zOffset - zMin) / zInc);
						int k = (int)kF;
						if(kF - k > 0.5f)
						{
							k++;
						}
						if(k != currentK)
						{
							planeData = (byte[])readPlaneData(ZDIR, getZCoor(k));
							currentK = k;
						}
						if(planeData == null)
						{
							value = 0;
						}
						else
						{
							byte signedByteValue = planeData[row * nCols + col];
							value = StsMath.signedByteToUnsignedInt(signedByteValue);
						}
						colors[surfaceRow][surfaceCol] = colorscaleColors[value];
					}
				}
			}
		}
		timer.stopPrint("get2dColorArray");
		return colors;
	}

	public byte[] getByteArray(StsRotatedGridBoundingBox surfaceBoundingBox, float[][] z, float zOffset)
	{
		int nSurfaceRows = surfaceBoundingBox.nRows;
		int nSurfaceCols = surfaceBoundingBox.nCols;

		byte[] bytes = new byte[nSurfaceRows * nSurfaceCols];

		int value;
		int currentK = -1;
		byte[] planeData = null;

		if(timer == null)
		{
			timer = new StsTimer();
		}
		timer.start();

		int surfaceRowMin = surfaceBoundingBox.getNearestBoundedRowCoor(yMin);
		int surfaceRowMax = surfaceBoundingBox.getNearestBoundedRowCoor(yMax);
		int surfaceColMin = surfaceBoundingBox.getNearestBoundedColCoor(xMin);
		int surfaceColMax = surfaceBoundingBox.getNearestBoundedColCoor(xMax);

		int n = 0;
		for(int surfaceRow = 0; surfaceRow < nSurfaceRows; surfaceRow++)
		{
			if(surfaceRow < surfaceRowMin || surfaceRow > surfaceRowMax)
			{
				continue;
			}

			int row = surfaceRow - surfaceRowMin;
			for(int surfaceCol = 0; surfaceCol < nSurfaceCols; surfaceCol++)
			{
				if(surfaceCol < surfaceColMin || surfaceCol > surfaceColMax)
				{
					continue;
				}

				int col = surfaceCol - surfaceColMin;

				float kF = (float)((z[surfaceRow][surfaceCol] + zOffset - zMin) / zInc);
				int k = (int)kF;
				if(kF - k > 0.5f)
				{
					k++;
				}
				if(k != currentK)
				{
					planeData = (byte[])readPlaneData(ZDIR, getZCoor(k));
					currentK = k;
				}
				if(planeData != null)
				{
					bytes[surfaceRow * nSurfaceCols + surfaceCol] = planeData[row * nCols + col];
				}
			}
		}
		timer.stopPrint("getByteArray");
		return bytes;
	}

	/** Display cursor is organized in rows and columns.  For an X-plane, display rows
	 *  are across the Y-axis and display columns are up the Z-axis so rows are from yMin to yMax
	 *  and columns are from zMin to zMax, for example.
	 */
	public float[][] getCursorDisplayRange(int dir)
	{
		switch(dir)
		{
			case XDIR:
				return new float[][]
					{
					{
					yMin, yMax}
					,
					{
					zMin, zMax}
				};
			case YDIR:
				return new float[][]
					{
					{
					xMin, xMax}
					,
					{
					zMin, zMax}
				};
			case ZDIR:
				return new float[][]
					{
					{
					xMin, xMax}
					,
					{
					yMin, yMax}
				};
			default:
				return null;
		}
	}

	/** Data on cursor is organized in rows and columns.  For an X-plane, data
	 *  increases down the Z axis and across the Y-axis so rows are from zMin to zMax
	 *  and columns are from yMin to yMax, for example.
	 */
	public float[][] getCursorDataRange(int dir)
	{
		switch(dir)
		{
			case XDIR:
				return new float[][]
					{
					{zMin, zMax},
					{yMin, yMax}
				};
			case YDIR:
				return new float[][]
					{
					{zMin, zMax},
					{xMin, xMax}
				};
			case ZDIR:
				return new float[][]
					{
					{xMin, xMax},
					{yMin, yMax}
				};
			default:
				return null;
		}
	}

	/** Map 2d coordinates isVisible on a cursor section to coordinates in which
	 *  data is actually organized.
	 */
	static public int[] getCursor2dCoorDataIndexes(int dir, boolean axesFlipped)
	{
		switch(dir)
		{
			case XDIR:
				return new int[]
					{1, 0};
			case YDIR:
				return new int[]
					{1, 0};
			case ZDIR:
				if(!axesFlipped)
				{
					return new int[]
						{0, 1};
				}
				else
				{
					return new int[]
						{1, 0};
				}
			default:
				return null;
		}
	}

	/** Given that point coordinates are X, Y, and Z, return Z,Y for X-plane,
	 *  Z,X for Y-plane, and X,Y for Z-plane.
	 */
	public int[] getCursorCoorDisplayIndexes(int dir, boolean axesFlipped)
	{
		switch(dir)
		{
			case XDIR:
				return new int[]
					{
					YDIR, ZDIR};
			case YDIR:
				return new int[]
					{
					XDIR, ZDIR};
			case ZDIR:

				if(!axesFlipped)
				{
					return new int[]
						{
						XDIR, YDIR};
				}
				else
				{
					return new int[]
						{
						YDIR, XDIR};
				}
			default:
				return null;
		}
	}

	public boolean canExport()
	{
		return true;
	}

	public boolean export()
	{
		try
		{
			if(isArchived)
			{
				StsDirectorySelectionDialog exportDialog = new StsDirectorySelectionDialog((Frame)currentModel.win3d,
					"Export Archived VSP File", segyFilename, true);
				exportDialog.setVisible(true);
				if(!exportDialog.getDoProcess())
					return false;
                exportDirectory = exportDialog.getCurrentDirectory();
				progressBarDialog = StsProgressBarDialog.constructor(currentModel.win3d, "Export Ancillary Data", false);
				progressBarDialog.setLabelText("Exporting " + getName());
			}
			else
			{
                new StsMessage(currentModel.win3d, StsMessage.WARNING, "Archived version not available, will re-create from processed data.");
				StsSeismicExportPanel.createDialog(currentModel, this, "VSP Export", false);
				return true;
			}
		}
		catch(Exception e)
		{
			new StsMessage(currentModel.win3d, StsMessage.ERROR, "Problem exporting VSP data");
		}

		Runnable runExport = new Runnable()
		{
			public void run()
			{
				exportData();
			}
		};
		Thread exportThread = new Thread(runExport);
		exportThread.start();
		return true;
	}

	private void exportData()
	{
		String outputName = null;
		try
		{
			String inputName = currentModel.getProject().getRootDirString() +
				currentModel.getProject().getArchiveDirString() + segyFilename;
			outputName = exportDirectory + File.separator + segyFilename;
			StsFile.copy(inputName, outputName);
		}
		catch(Exception e)
		{
			new StsMessage(currentModel.win3d, StsMessage.ERROR, "Error exporting " + outputName);
		}
	}

	/** Return number of rows and cols on isVisible 2D plane from 3D cube of nRows, nCols, nCroppedSlices.
	 *  Return nCroppedSlices,nRows for X-plane, nCroppedSlices,nCols for Y-plane, and nRows,nCols for Z-plane.
	 */
	public int[] getCursorDisplayNRowCols(int dir)
	{
		switch(dir)
		{
			case XDIR:
				return new int[]
					{
					nSlices, nRows};
			case YDIR:
				return new int[]
					{
					nSlices, nCols};
			case ZDIR:
				return new int[]
					{
					nRows, nCols};
			default:
				return null;
		}
	}

	/** Return number of rows and cols on data plane from 3D cube of nRows, nCols, nCroppedSlices.
	 *  Texture goes down first trace on vertical planes and across first row on horizontal plane.
	 *  Return nCols, nCroppedSlices for X-plane, nRows,nCroppedSlices for Y-plane, and nRows,nCols for Z-plane.
	 */
	public int[] getCursorDataNRowCols(int dir)
	{
		switch(dir)
		{
			case XDIR:
				return new int[]
					{
					nRows, nSlices};
			case YDIR:
				return new int[]
					{
					nCols, nSlices};
			case ZDIR:
				return new int[]
					{
					nRows, nCols};
			default:
				return null;
		}
	}

	public Object readPlaneData(int dir, float dirCoordinate)
	{
		try
		{
			return readPlaneData(1, 0);
		}
		catch(Exception e)
		{
			StsException.outputException("StsVsp.readPlaneData() failed.", e, StsException.WARNING);
			return null;
		}
	}

	public float getAmplitudeAt(int nTrace, float time)
	{
		nTrace = StsMath.minMax(nTrace, 0, nCols - 1);
		int nSlice = this.getNearestBoundedSliceCoor(time);
		return readValue(nTrace, nSlice);
	}

	/*	public byte[] getByteData()
	 {
	  if(byteData != null)return byteData;
	  getScaledFloatData();
	  return byteData;
	 }

	 public byte[] getByteData(boolean flipped)
	 {
	  if(!flipped)
	   return getByteData();
	  else
	  {
	   byte[] d = getByteData();
	   byte[] rdata = new byte[d.length];
	   for(int i = 0; i < d.length; i++)
		rdata[i] = (byte) - byteData[i];
	   return rdata;
	  }
	 }
	 */

	public float[] getFloatData()
	{
		if(scaledFloatData != null)return scaledFloatData;
		floatData = readFloatData();
		return floatData;
	}

	/*	public float[] getFloatData(boolean flipped)
	 {
	  if(!flipped)
	   return getFloatData();
	  else
	  {
	   float[] d = getFloatData();
	   float[] rdata = new float[d.length];
	   for(int i = 0; i < d.length; i++)
		rdata[i] = (float) - byteData[i];
	   return rdata;
	  }
	 }
	 */
	public float[] getScaledFloatData()
	{
		if(scaledFloatData != null)return scaledFloatData;
		floatData = getFloatData();
		if(floatData == null)return null;
		scaledFloatData = StsSeismicFilter.filter(filterProperties.PRESTACK, floatData, nCols, nSlices, getZInc(), filterProperties, agcProperties, dataMin, dataMax);
		StsSeismicFilter.normalizeAmplitude(scaledFloatData, nCols*nSlices);
		return scaledFloatData;
	}

    public byte[] getData()
    {
        if(byteData != null) return byteData;
        getScaledFloatData();
        if(scaledFloatData == null) return null;
        byteData = StsMath.floatsToUnsignedBytes254(scaledFloatData, dataMin, dataMax);
        return byteData;
    }

    public byte[] readData()
	{
		try
		{
			DataInputStream dis = StsFile.constructDIS(stsDirectory + rowCubeFilename);
			int nValues = nCols * nSlices;
			byte[] values = new byte[nValues];
			dis.read(values);
			dis.close();
			return values;
		}
		catch(Exception e)
		{
			StsException.systemError("StsVsp.readData() failed to find file " + stsDirectory + rowCubeFilename);
			return null;
		}
	}

	public float[] readFloatData()
	{
		try
		{
			StsMappedFloatBuffer floatBuffer = StsMappedFloatBuffer.openRead(stsDirectory, rowFloatFilename);
			int nValues = nCols * nSlices;
			if(!floatBuffer.map(0, nValues)) return null;
			float[] values = new float[nValues];
			floatBuffer.get(values);
			floatBuffer.close();
			return values;
		}
		catch(Exception e)
		{
			StsException.systemError("StsVsp.readFloatData() failed to find file " + stsDirectory + rowFloatFilename);
			return null;
		}
	}

	public float readValue(int nTrace, int nSlice)
	{
		float ret = StsParameters.nullValue;
		if(nTrace < 0 || nTrace >= nCols)return StsParameters.nullValue;
		if(nSlice < 0 || nSlice >= nSlices)return StsParameters.nullValue;
		float[] f = getScaledFloatData();
		if(f != null)
			ret = f[nTrace * nSlices + nSlice];
		return ret;
	}

	public float readRawValue(int nTrace, int nSlice)
	{
		try
		{
			StsMappedFloatBuffer floatBuffer = StsMappedFloatBuffer.openRead(stsDirectory, rowFloatFilename);
			if(nTrace < 0 || nTrace >= nCols)return StsParameters.nullValue;
			if(nSlice < 0 || nSlice >= nSlices)return StsParameters.nullValue;
			int nValues = nCols * nSlices;
			if(!floatBuffer.map(0, nValues)) return StsParameters.nullValue;
			floatBuffer.position(nTrace * nSlices + nSlice);
			float value = floatBuffer.getFloat();
			floatBuffer.close();
			return value;
		}
		catch(Exception e)
		{
			StsException.systemError("StsVsp.readFloatData() failed to find file " + stsDirectory + rowFloatFilename);
			return StsParameters.nullValue;
		}
	}

	public boolean checkDataFile() throws FileNotFoundException
	{
		try
		{
			File file = new File(stsDirectory + rowCubeFilename);
			return file.exists();
		}
		catch(Exception e)
		{
			return false;
		}
	}

	public float getScaledValue(byte byteValue)
	{
		float f = (float)StsMath.signedByteToUnsignedInt(byteValue);
		return dataMin + (f / 254) * (dataMax - dataMin);
	}

	public byte getByteValueFromFloat(float value)
	{
		return StsMath.unsignedIntToUnsignedByte((int)(254 * (value - dataMin) / (dataMax - dataMin)));
	}

	public final boolean isByteValueNull(byte byteValue)
	{
		return byteValue == -1;
	}

	/** in byte range 0 to 254, the zero crossing is at this value */
	public float getUnsignedByteAverage()
	{
		return 127;
	}

	public synchronized void addItemListener(ItemListener l)
	{
		Vector v = itemListeners == null ? new Vector(2) : (Vector)itemListeners.clone();
		if(!v.contains(l))
		{
			v.addElement(l);
			itemListeners = v;
		}
	}

	public synchronized void removeItemListener(ItemListener l)
	{
		if(itemListeners != null && itemListeners.contains(l))
		{
			Vector v = (Vector)itemListeners.clone();
			v.removeElement(l);
			itemListeners = v;
		}
	}

	protected void fireItemStateChanged(ItemEvent e)
	{
		if(itemListeners != null)
		{
			Vector listeners = itemListeners;
			int count = listeners.size();
			for(int i = 0; i < count; i++)
			{
				((ItemListener)listeners.elementAt(i)).itemStateChanged(e);
			}
		}
	}

	public void accumulateHistogram(int bindex)
	{
		if(bindex > 254)
		{
			bindex = 254;
		}
		if(bindex < 0)
		{
			bindex = 0;
		}
		dataCnt[bindex] = dataCnt[bindex] + 1;
		ttlHistogramSamples++;
	}

	public void calculateHistogram()
	{
		for(int i = 0; i < 255; i++)
		{
			dataHist[i] = ((float)dataCnt[i] / (float)ttlHistogramSamples) * 100.0f;
		}
	}

	public void clearHistogram()
	{
		for(int i = 0; i < 255; i++)
		{
			dataCnt[i] = 0;
			dataHist[i] = 0.0f;
		}
		ttlHistogramSamples = 0;
	}

	public void drawTextureTileSurface(StsTextureTile tile, GL gl, int dir, boolean is3d, int nTile)
	{
		byte projectZDomain = StsObject.getCurrentModel().getProject().getZDomain();
		byte volumeZDomain = getZDomain();

		if(projectZDomain == StsParameters.TD_TIME)
		{
			if(volumeZDomain == StsParameters.TD_TIME)
				drawTextureTileTimeSurface(tile, gl, is3d);
		}
		else if(projectZDomain == StsParameters.TD_DEPTH)
		{
			// volumeZDomain is TD_TIME
			StsModel model = StsObject.getCurrentModel();
			if(volumeZDomain == StsParameters.TD_DEPTH) // seismic already in depth, don't need to convert so draw as if in time
				drawTextureTileTimeSurface(tile, gl, is3d);
			else
			{
				StsSeismicVelocityModel velocityVolume = model.getProject().velocityModel;
				if(velocityVolume == null)return;
				drawTextureTileDepthSurface(velocityVolume, tile, gl, dir);
			}
		}
	}

	public void drawTextureTileTimeSurface(StsTextureTile tile, GL gl, boolean is3d)
	{
		if(is3d)
			tile.drawQuadSurface3d(gl, StsProject.TD_TIME);
		else
			tile.drawQuadSurface2d(gl);
	}

	public void drawTextureTileDepthSurface(StsSeismicVelocityModel velocityModel, StsTextureTile tile, GL gl, int dir)
	{
		float cursorXInc, cursorYInc;
		if(dir == StsCursor3d.ZDIR)
		{
			return;
		}

		gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
		double rowTexCoor = tile.minRowTexCoor;
		double dRowTexCoor = tile.dRowTexCoor;
		double dColTexCoor = tile.dColTexCoor;
		double[] xyz = tile.xyzPlane[0];
		double x1 = xyz[0];
		double y1 = xyz[1];
		double t1 = xyz[2];
		int volumeRow = this.getNearestBoundedRowCoor((float)y1);
		int volumeCol = getNearestBoundedColCoor((float)x1);
		StsSeismicVolume velocityVolume = velocityModel.getVelocityVolume();
		if(velocityVolume == null)return;
		float depthMin = velocityModel.depthDatum;
		int volumeRowInc = 0;
		int volumeColInc = 0;
		if(dir == StsCursor3d.XDIR)
		{
			cursorXInc = 0;
			cursorYInc = getYInc();
			volumeRowInc = 1;
		}
		else // dirNo == StsCursor3d.YDIR
		{
			cursorXInc = getXInc();
			cursorYInc = 0;
			volumeColInc = 1;
		}
		double tInc = getZInc();
		for(int row = tile.croppedRowMin + 1; row <= tile.croppedRowMax; row++, rowTexCoor += dRowTexCoor)
		{
			double x0 = x1;
			double y0 = y1;
			x1 += cursorXInc;
			y1 += cursorYInc;

			gl.glBegin(GL.GL_QUAD_STRIP);

			double colTexCoor = tile.minColTexCoor;
			double t = t1 + tile.croppedColMin * tInc;

			for(int col = tile.croppedColMin; col <= tile.croppedColMax; col++, t += tInc, colTexCoor += dColTexCoor)
			{
				float v0 = velocityVolume.getValue(volumeRow, volumeCol, col);
				float z0 = (float)(v0 * t + depthMin);
				gl.glTexCoord2d(colTexCoor, rowTexCoor);
				gl.glVertex3d(x0, y0, z0);
				float v1 = velocityVolume.getValue(volumeRow + volumeRowInc, volumeCol + volumeColInc, col);
				float z1 = (float)(v1 * t + depthMin);
				gl.glTexCoord2d(colTexCoor, rowTexCoor + dRowTexCoor);
				gl.glVertex3d(x1, y1, z1);
			}
			gl.glEnd();
			volumeRow += volumeRowInc;
			volumeCol += volumeColInc;
		}
	}

	public void deleteDisplayList(GL gl)
	{
		if(displayListNum > 0)
		{
			gl.glDeleteLists(displayListNum, 1);
			displayListNum = 0;
		}
	}

	private void setColorListNum(int colorListNum)
	{
		this.colorListNum = colorListNum;
	}

	final public byte floatToSignedByte254(float value, float scale)
	{
		int i = Math.round(dataMin + (value - dataMin) * scale);
		if(i >= 255)i = 254;
		if(i < 0)i = 0;
		return(byte)i;
	}

	public StsWell getWell()
	{
		return well;
	}

	public void setWell(StsWell well)
	{
        this.well = well;
		fieldChanged("well", well);
	}

	public void objectPropertiesChanged()
	{
		scaledFloatData = null;
		textureChanged = true;

	}

	public void display(StsGLPanel glPanel)
	{
        StsGLPanel3d glPanel3d = (StsGLPanel3d)glPanel;
		GL gl = glPanel3d.getGL();

		if(well == null)return;
		if(!well.getIsVisible())return;
		/*
		 if(isChanged())
		 {
		  scaledFloatData = null;
//			byteData = null;
		  //if(vspView != null && vspView.getWellVspDisplayPanel() != null)
		  //	vspView.getWellVspDisplayPanel().setTextureChanged(true);
		  textureChanged = true;
		 }
		 */
		if(!displayVSPs)return;
		if(textureTiles == null)
		{
			if(!initializeTextureTiles(glPanel3d, gl))return;
			textureChanged = true;

		}
		StsView3d view3d = (StsView3d)glPanel3d.getView();
		azimuth = view3d.azimuth;

		float[] xyz = well.getXYZPointAtZorT(zMin, true).v;
		gl.glPushMatrix();
		gl.glTranslatef(xyz[0], xyz[1], xyz[2]);
		gl.glRotated( -azimuth, 0.0, 0.0, -1.0);

        glPanel3d.setViewShift(gl, -StsGraphicParameters.edgeShift);
        if(displayVAR)
			displayTexture(glPanel3d, gl);
        glPanel3d.resetViewShift(gl);
		if(getWiggleDisplayProperties().getDisplayWiggles())
			displayWiggleTraces(glPanel3d, gl);
		gl.glPopMatrix();
	}

	private void displayTexture(StsGLPanel3d glPanel3d, GL gl)
	{
		gl.glDisable(GL.GL_LIGHTING);
		gl.glEnable(GL.GL_TEXTURE_2D);
		gl.glEnable(GL.GL_BLEND);
		gl.glTexEnvf(GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_MODULATE);
		gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);
		gl.glShadeModel(GL.GL_FLAT);

		if(colorListChanged)textureChanged = true;

//		glPanel3d.setViewShift(gl, -StsGraphicParameters.edgeShift);

        if(textureTiles.shaderChanged()) textureChanged = true;

		if(isPixelMode != getIsPixelMode())
		{
			textureChanged = true;
			isPixelMode = !isPixelMode;
		}

	    setGLColorList(gl, false, textureTiles.shader);

		if(textureChanged)
		{
			textureTiles.displayTiles(this, gl, isPixelMode, getScaledFloatData(), nullByte);
			textureChanged = false;
			colorListChanged = false;
			// jbw
			//if(vspView != null)vspView.repaint();
		}
		else
			textureTiles.displayTiles(this, gl, isPixelMode, (byte[])null, nullByte);
//		glPanel3d.resetViewShift(gl);
		gl.glDisable(GL.GL_TEXTURE_2D);
		gl.glDisable(GL.GL_BLEND);
		gl.glEnable(GL.GL_LIGHTING);

		if(getTextureTiles().shader != StsJOGLShader.NONE) StsJOGLShader.disableARBShader(gl);
	}

    private boolean initializeTextureTiles(StsGLPanel3d glPanel3d, GL gl)
	{
		// axisRanges have origin at lower-left (xmin, zmax) and first axis is horizontal, second is vertical up.
		// textureRanges have origin at upper-left (xmin, zmin) and first axis is vertical down and second is horizontal across.
		panelRanges = new float[2][2];
		panelRanges[0][0] = zMin;
		panelRanges[0][1] = zMax;
		panelRanges[1][0] = -panelWidth;
		panelRanges[1][1] = panelWidth;
		if(textureTiles != null) this.deleteTexturesAndDisplayLists(gl);
		textureTiles = StsTextureTiles.constructor(currentModel, this, nCols, nSlices, true, panelRanges);
		if(textureTiles == null)return false;
		//data = getData(); // not yet! wait for shader set-up
		//textureTiles.addData(gl, data);
		//textureChanged = true;

		return true;
	}

	public void drawTextureTileSurface(StsTextureTile tile, GL gl, boolean is3d)
	{
		tile.drawQuadSurface3d(gl, StsProject.TD_TIME);
	}

    public void deleteTexturesAndDisplayLists(GL gl)
    {
       if(textureTiles == null)return;
       textureTiles.deleteTextures(gl);
       textureChanged = true;
    }

	public boolean textureChanged()
	{
        textureChanged = true;
        return true;
    }

    public boolean dataChanged()
    {
        return textureChanged();
    }
    
    public void geometryChanged()
    {
    }
    
    private void displayWiggleTraces(StsGLPanel3d glPanel3d, GL gl)
	{
		float[] data = getScaledFloatData();
		if(data == null)return;

		StsPoint wellTopPoint = well.getXYZPointAtZorT(zMin, true);
		StsPoint downPoint = new StsPoint(wellTopPoint);
		downPoint.setZ(zMax);
		StsPoint nearestPoint = currentModel.getGlPanel3d().getPointOnLineNearestViewLine(wellTopPoint, downPoint);
		if(wiggleDebug)System.out.println("displayWiggles Nearest point: " + nearestPoint.toString());
		// Never display wiggles in 3D unless close enough to make sense
		double[] wellScreenPoint = currentModel.getGlPanel3d().getScreenCoordinates(nearestPoint);
		wellScreenPoint[0] += 100;
        //TODO instead, use getWorldCoordinates which returns double[] 
        StsPoint sidewaysPoint = currentModel.getGlPanel3d().getWorldCoordinatesPoint(wellScreenPoint);
		if(wiggleDebug)System.out.println("displayWiggles sidewaysPoint (100 pixels away in screen x): " + sidewaysPoint.toString());
		float lengthPerPixel = nearestPoint.distance(sidewaysPoint) / 100;
		double panelPixels = 2 * panelWidth / lengthPerPixel;
		if(wiggleDebug)
		{
			System.out.println("lengthPerPixel " + lengthPerPixel + " PanelPixels " + panelPixels);
			System.out.println("pixelsPerTrace " + panelPixels / nCols + " wiggleToPixelRatio " + wiggleToPixelRatio);
		}
		if(panelPixels / nCols < wiggleToPixelRatio)return;

		float traceWidth = (float)(2 * panelWidth) / nCols;
//		System.out.println("traceWidth " + traceWidth);
		//float horizScale = 3 * traceWidth / 254;
		float horizScale = traceWidth / 2;
		if(getWiggleDisplayProperties().getWiggleReversePolarity())
			horizScale = -horizScale;
		gl.glDisable(GL.GL_LIGHTING);
		gl.glLineWidth(0.5f);
		glPanel3d.setViewShift(gl, StsGraphicParameters.edgeShift);
		gl.glDisable(GL.GL_LIGHTING);
		int n = 0;
		float x = -panelWidth + traceWidth / 2;
		for(int t = 0; t < nCols; t++, x += traceWidth, n += nSlices)
		{
			float z = zMin;
			drawFilledWiggleTraces(gl, data, horizScale, n, x);
			drawWiggleTraces(gl, data, horizScale, n, x);
		}
		gl.glEnable(GL.GL_LIGHTING);
		glPanel3d.resetViewShift(gl);
	}

	private void drawFilledWiggleTraces(GL gl, float[] data, float horizScale, int n, float x0)
	{

		try
		{
			if(!getWiggleDisplayProperties().hasFill())return;
			StsColor plusColor = getWiggleDisplayProperties().getWigglePlusColor();
			StsColor minusColor = getWiggleDisplayProperties().getWiggleMinusColor();
			float t1 = zMin;
			float a1 = horizScale * data[n];
			boolean b1 = a1 >= 0;
			if(b1)
                plusColor.setGLColor(gl);
			else
				minusColor.setGLColor(gl);

			gl.glBegin(GL.GL_QUAD_STRIP);
			gl.glVertex3f(x0, 0.0f, t1);
			gl.glVertex3f(x0 + a1, 0.0f, t1);

			for(int s = n + 1; s < n + nSlices; s++)
			{
				float t0 = t1;
				t1 = t1 + zInc;
				float a0 = a1;
				a1 = horizScale * data[s];
				boolean b0 = b1;
				b1 = a1 >= 0;
				if(b0 && b1 || !b0 && !b1)
				{
					gl.glVertex3f(x0, 0.0f, t1);
					gl.glVertex3f(x0 + a1, 0.0f, t1);
				}
				else
				{
					float tm = t0 + a0 * (t1 - t0) / (a0 - a1);
					gl.glVertex3f(x0, 0.0f, tm);
					gl.glVertex3f(x0, 0.0f, tm);
                    if(b1)
                        plusColor.setGLColor(gl);
                    else
                        minusColor.setGLColor(gl);
					gl.glVertex3f(x0, 0.0f, tm);
					gl.glVertex3f(x0, 0.0f, tm);
					gl.glVertex3f(x0, 0.0f, t1);
					gl.glVertex3f(x0 + a1, 0.0f, t1);
				}
			}
		}
		catch(Exception e)
		{
			StsException.outputException("StsVSP.displayFilledWiggleTraces() failed.", e, StsException.WARNING);
		}
		finally
		{
			gl.glEnd();
		}
	}

	private void drawWiggleTraces(GL gl, float[] data, float horizScale, int n, float x0)
	{
		try
		{
			if(displayVAR)
				StsColor.BLACK.setGLColor(gl);
			else
				StsColor.WHITE.setGLColor(gl);
			float z = zMin;
			gl.glBegin(GL.GL_LINE_STRIP);
			for(int s = n; s < n + nSlices; s++, z += zInc)
				gl.glVertex3d(x0 + horizScale * (data[s]), 0.0, z);
		}
		catch(Exception e)
		{
			StsException.outputException("StsVSP.displayWiggleTraces() failed.", e, StsException.WARNING);
		}
		finally
		{
			gl.glEnd();
		}
	}

	public boolean isChanged()
	{
		if(getWiggleDisplayProperties().isChanged())return true;
		if(filterProperties.isChanged())return true;
		if(agcProperties.isChanged()) return true;
        if(datumProperties.isChanged()) return true;
		return false;
	}

	public void displayVspPropertiesDialog()
	{
		displayVspPropertiesDialog(currentModel.win3d);
	}

	public void displayVspPropertiesDialog(Frame frame)
	{
		new StsOkApplyCancelDialog(frame, new StsDialogFace[]
								   {getWiggleDisplayProperties(), agcProperties, datumProperties, filterProperties}, getName() + " Vsp Properties", false);
//		currentModel.viewObjectRepaint(this);
	}

    public double[] getAttributeArrayInMdepth(String attributeName, byte domain) throws FileNotFoundException
    {
        attributeFilename = attributeFilePrefix + stemname;
        attributeDomain = domain;
        attributeInMdepth = true;
        double[] atts = getAttributeArray(stsDirectory, attributeFilename, attributeName);
        attributeInMdepth = false;
        return atts;
	}

	public double[] getAttributeArray(String attributeName) throws FileNotFoundException
	{
		attributeFilename = attributeFilePrefix + stemname;
		return getAttributeArray(stsDirectory, attributeFilename, attributeName);
	}

	private double[] getAttributeArray(String directoryname, String filename, String attributeName) throws FileNotFoundException
	{
		RandomAccessFile raf = null;
		int attributeIndex = this.getAttributeIndex(attributeName);
		if(attributeIndex == -1)
		{
			//new StsMessage(currentModel.win3d, StsMessage.WARNING, "Failed to find " + attributeName + " attribute in attribute names list.");
			return null;
		}
		try
		{
			raf = new RandomAccessFile(directoryname + filename, "r");
			long offset = (long)attributeIndex * 8 * nCols;
			raf.seek(offset);
			byte[] bytes = new byte[8 * nCols];
			raf.readFully(bytes);
			ByteArrayInputStream bais = new ByteArrayInputStream(bytes);
			DataInputStream dis = new DataInputStream(bais);
			double[] attributes = new double[nCols];
			for(int n = 0; n < nCols; n++)
			{
                attributes[n] = dis.readDouble();
                if((well != null) && (attributeInMdepth))
                {
                    if(attributeDomain == StsParameters.TD_TIME)
                        attributes[n] = (double) well.getMDepthFromTime((float) attributes[n]);
                    else
                        attributes[n] = (double) well.getMDepthFromDepth((float) attributes[n]);
                }
			}
			return attributes;
		}
		catch(Exception e)
		{
			StsException.systemError("StsVsp.getAttributeArray failed.");
			return null;
		}
		finally
		{
			try
			{
				raf.close();
			}
			catch(Exception e)
			{}
		}
	}

    public void setBornDate(String born)
    {
        if(!StsDateFieldBean.validateDateInput(born))
        {
            bornField.setValue(StsDateFieldBean.convertToString(bornDate));
            return;
        }
        super.setBornDate(born);
    }
    public void setDeathDate(String death)
    {
        if(!StsDateFieldBean.validateDateInput(death))
        {
            deathField.setValue(StsDateFieldBean.convertToString(deathDate));
            return;
        }
        super.setDeathDate(death);
    }
    // Demo Purposes
//    public String getOwner() { return "Pemex"; }
//    public String getAcquisitionDate() { return "Jan. 12, 2003"; }
//    public String getField() { return "Cantrell Cplx"; }
//    public String getAcquiredBy() { return "Schlumberger"; }
//    public String getProcessor() { return "Jaquar Expl."; }
//    public String getProcessDate() { return "Jan. 31, 2006"; }
	/*
	 double[] topScreenCoors = glPanel3d.getScreenCoordinates(topPoint.v);
	 topScreenCoors[1] -= 1;
	 StsPoint topLeftPoint = glPanel3d.getWorldCoordinatesPoint(topScreenCoors);
	 StsPoint vector = new StsPoint(3);
	 vector.subPoints(topLeftPoint, topPoint);
	 float length = vector.length();
	 float scaleFactor = 100/length;
	 */

    public StsWiggleDisplayProperties getWiggleDisplayProperties()
    {
        return vspClass.getWiggleDisplayProperties();
    }

    public Class getDisplayableClass() { return StsVsp.class; }

    public StsTextureTiles getTextureTiles()
    {
        return textureTiles;
    }
}