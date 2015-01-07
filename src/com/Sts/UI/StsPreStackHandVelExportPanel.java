package com.Sts.UI;

import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.Types.PreStack.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.util.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

public abstract class StsPreStackHandVelExportPanel extends StsJPanel implements StsDialogFace
{
	protected String exportHeader = null;
	public byte exportDomainType = StsParameters.TD_TIME;
	protected float rowNumMin, rowNumMax, rowNumInc;
	protected float colNumMin, colNumMax, colNumInc;
	protected float zMin, zMax, zInc;
	protected String velocityUnits = StsParameters.VEL_FT_PER_MSEC;
	public String allOrSomeString = ALL_PROFILES;
	protected float unitsScale = 1.0f;
	protected byte nullType = NULL_ZERO;
	protected byte velocityType = StsParameters.SAMPLE_TYPE_VEL_RMS;
	public String textHeader = null;

	static public final String ALL_PROFILES = "All Hand Picked";
	static public final String SOME_PROFILES = "Some Hand Picked";
	static public final String INTERPOLATED_PROFILES = "Interpolated Profiles";

	protected transient StsPreStackVelocityModel volume = null;
	public transient StsPreStackLineSet lineSet = null;
	protected transient JLabel titleLabel = new JLabel("Velocity Profile Export");
	protected transient String[] profileStrings = new String[]
													{ALL_PROFILES, SOME_PROFILES, INTERPOLATED_PROFILES};
	protected transient StsButtonListFieldBean profileButtons = new StsButtonListFieldBean(this, "allOrSomeProfiles", null, profileStrings, true);
	protected transient StsGroupBox exportParametersGroupBox = new StsGroupBox("File Export Parameters");
	protected transient StsStringFieldBean exportNameBean;
	protected transient StsComboBoxFieldBean nullValueBean;
	protected transient StsGroupBox subVolumeGroupBox = new StsGroupBox("Specify Area to Export");
	protected transient JLabel minLabel = new JLabel("Minimum");
	protected transient JLabel maxLabel = new JLabel("Maximum");
	protected transient JLabel intervalLabel = new JLabel("Increment");
	protected transient StsFloatFieldBean rowNumMinBean;
	protected transient StsFloatFieldBean rowNumMaxBean;
	protected transient StsFloatFieldBean rowNumIncBean;
	protected transient StsFloatFieldBean colNumMinBean;
	protected transient StsFloatFieldBean colNumMaxBean;
	protected transient StsFloatFieldBean colNumIncBean;
	protected transient StsFloatFieldBean zMinBean;
	protected transient StsFloatFieldBean zMaxBean;
	protected transient StsFloatFieldBean zIncBean;
	protected transient StsGroupBox unitsGroupBox = null;
	transient static String group = "handVelExport";
	transient static String format = "txt";
	protected transient String exportName = null;
	
	public String[] pathnames = null; /** list of HandVel filenames to export */

	// these members are associated with the unitsGroupBox which is added to the seismicExportPanel
	transient StsFloatFieldBean unitsScaleBean;

	protected transient StsComboBoxFieldBean velocityTypeBean;
	protected transient StsComboBoxFieldBean exportDomainBean;
	protected transient StsParameterFileFieldBean templateBean;

	public final static String[] VEL_STRINGS = StsParameters.VEL_STRINGS;

	public final static byte NULL_NONE = 0;
	public final static byte NULL_MIN = 1;
	public final static byte NULL_MAX = 2;
	public final static byte NULL_ZERO = 3;

	public final static String NULL_NONE_STRING = "None";
	public final static String NULL_MIN_STRING = "Min";
	public final static String NULL_MAX_STRING = "Max";
	public final static String NULL_ZERO_STRING = "Zero";
	public final static String[] nullValueStrings = new String[]
																	{NULL_NONE_STRING, NULL_MIN_STRING, NULL_MAX_STRING, NULL_ZERO_STRING};

//	StsFloatFieldBean startInlineBean = new StsFloatFieldBean();
//	StsFloatFieldBean endInlineBean = new StsFloatFieldBean();
//	StsFloatFieldBean startCrosslineBean = new StsFloatFieldBean();
//	StsFloatFieldBean endCrosslineBean = new StsFloatFieldBean();

	protected transient StsButton textHdrEditBtn = new StsButton("Edit File Header", "Edit the default text header info.", this, "editFileHeader");
//	StsButton processBtn = new StsButton("Process", "Run the hand vel export process.", this, "process");

	protected transient StsModel model;

	transient public int okCancelType;

	public final static byte PROCESS = 0;
	public final static byte CANCELED = 1;

	protected transient boolean canceled = false;

	public StsPreStackHandVelExportPanel(StsModel model, StsPreStackVelocityModel volume)
	{
		this.volume = volume;
		this.lineSet = volume.getSeismicVolume();
		this.velocityUnits = model.getProject().getVelocityUnits();
		this.model = model;
		exportName = StsStringUtils.cleanString(volume.getName());
		this.unitsGroupBox = unitsGroupBox;
		try
		{
			buildBeans();
			assembleGroupBoxes();
			assemblePanel();
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
	}

	protected void buildBeans()
	{
		templateBean = new StsParameterFileFieldBean(this, group, format);
		templateBean.setComboBoxToolTipText("Select existing handvel export template");
		templateBean.setNewButtonToolTipText("Create a new handvel export template and save it");

		exportNameBean = new StsStringFieldBean(this, "exportName", "Export name:");
		exportNameBean.setToolTipText("Exports to file directory/name.handvel");
		velocityTypeBean = new StsComboBoxFieldBean(this, "velocityTypeString", "Velocity type:", StsParameters.VEL_STRINGS);
		exportDomainBean = new StsComboBoxFieldBean(this, "exportDomainString", "Export type:", StsParameters.TD_STRINGS);
		exportDomainBean.setEditable(false);
        nullValueBean = new StsComboBoxFieldBean(this, "nullValueString", "Null value:", this.nullValueStrings);

		rowNumMin = volume.getRowNumMin();
		rowNumMax = volume.getRowNumMax();
		rowNumInc = volume.getRowNumInc();
		colNumMin = volume.getColNumMin();
		colNumMax = volume.getColNumMax();
		colNumInc = volume.getColNumInc();
		zMin = volume.getZMin();
		zMax = volume.getZMax();
		zInc = volume.getZInc();
		rowNumMinBean = new StsFloatFieldBean(this, "rowNumMin", "Inline:");
		rowNumMaxBean = new StsFloatFieldBean(this, "rowNumMax");
		rowNumIncBean = new StsFloatFieldBean(this, "rowNumInc");
		colNumMinBean = new StsFloatFieldBean(this, "colNumMin", "Crossline:");
		colNumMaxBean = new StsFloatFieldBean(this, "colNumMax");
		colNumIncBean = new StsFloatFieldBean(this, "colNumInc");
		zMinBean = new StsFloatFieldBean(this, "zMin", "Time/Depth:");
		zMaxBean = new StsFloatFieldBean(this, "zMax");
		zIncBean = new StsFloatFieldBean(this, "zInc");

		setAllOrSomeProfiles(allOrSomeString);

//        vUnitsBean = new StsComboBoxFieldBean(this, "velocityUnits", "Velocity Units:", StsParameters.VEL_UNITS );
//        vUnitsBean.setToolTipText("What units would you like for the velocity output?");
	}

	protected abstract void assembleGroupBoxes() throws Exception;

	protected void assemblePanel()
	{
		titleLabel.setFont(new java.awt.Font("Serif", 1, 14));
		gbc.fill = GridBagConstraints.HORIZONTAL;
		add(titleLabel);
		add(profileButtons.getButtonPanel());
		add(exportParametersGroupBox);
		add(subVolumeGroupBox);
		add(unitsGroupBox);
		gbc.fill = GridBagConstraints.NONE;
		add(textHdrEditBtn);
		gbc.fill = GridBagConstraints.HORIZONTAL;
		//add(progressPanel);
	}

	private void setAllOrSomeProfiles(Object profileType)
	{
		allOrSomeString = (String)profileType;
		if (profileType == SOME_PROFILES)
		{
			rowNumIncBean.setEditable(false);
			rowNumMinBean.setEditable(true);
			rowNumMaxBean.setEditable(true);
			colNumIncBean.setEditable(false);
			colNumMinBean.setEditable(true);
			colNumMaxBean.setEditable(true);
			zIncBean.setEditable(false);
			zMinBean.setEditable(false);
			zMaxBean.setEditable(false);
		}
		else if (profileType == INTERPOLATED_PROFILES)
		{
			rowNumIncBean.setEditable(true);
			rowNumMinBean.setEditable(true);
			rowNumMaxBean.setEditable(true);
			colNumIncBean.setEditable(true);
			colNumMinBean.setEditable(true);
			colNumMaxBean.setEditable(true);
			zIncBean.setEditable(true);
			zMinBean.setEditable(true);
			zMaxBean.setEditable(true);
		}
		else if (profileType == ALL_PROFILES)
		{
			rowNumIncBean.setEditable(false);
			rowNumMinBean.setEditable(false);
			rowNumMaxBean.setEditable(false);
			colNumIncBean.setEditable(false);
			colNumMinBean.setEditable(false);
			colNumMaxBean.setEditable(false);
			zIncBean.setEditable(false);
			zMinBean.setEditable(false);
			zMaxBean.setEditable(false);
		}
	}

	protected Object getAllOrSomeProfiles()
	{
		return this.allOrSomeString;
	}

	public String getProfileType()
	{
		return allOrSomeString;
	}

	public void editFileHeader()
	{
		if (textHeader == null)
			textHeader = getFileHeader();

		StsTextAreaDialog volumeDialog = new StsTextAreaDialog(model.win3d, "HandVel Export Header", textHeader, 20, 80);
		volumeDialog.setVisible(true);
		textHeader = volumeDialog.getText();
	}

	public void dialogSelectionType(int type)
	{
		if (type == StsDialogFace.PROCESS || type == StsDialogFace.OK)
		{
			run();
			return;
		}
		if (type == StsDialogFace.CANCEL)
		{
			canceled = true;
		}
	}

	public String getVelocityUnitsString()
	{
		return this.velocityUnits;
	}

	public void setVelocityUnitsString(String unitsString)
	{
		this.velocityUnits = unitsString;
		unitsScale = 1.0f / model.getProject().calculateVelScaleMultiplier(unitsString);
		unitsScaleBean.setValue(unitsScale);
	}

    public float getUnitsScale() { return unitsScale; }

    public boolean run()
	{
		final StsDialogFace face = this;
		Runnable runExport = new Runnable()
		{
			public void run()
			{
				StsProgressCancelDialog dialog = null;
				try
				{
					dialog = new StsProgressCancelDialog(model.win3d, face, "Export Hand Vels.", "Export Progress", false);
					StsProgressPanel progressPanel = dialog.getProgressPanel();
					dialog.setVisible(true);
					exportHandVels(progressPanel);
				}
				finally
				{
					dialog.setProcessCompleted();
				}
			}
		};

		Thread exportThread = new Thread(runExport);
		exportThread.start();
		return true;
	}

	public abstract boolean exportHandVels(StsProgressPanel progressPanel);

	private boolean askedBeforeDontAskAgain = false;

	public StsAsciiFile setupExport(String filename, StsPreStackLine line) throws IOException
	{
		File file = new File(filename);
		if (file.exists() && ! askedBeforeDontAskAgain)
		{
			boolean overWrite = StsYesNoDialog.questionValue(model.win3d, "File " + filename + " already exists. Do you wish to overwrite it?");
			if (!overWrite)
				return null;
			askedBeforeDontAskAgain = true;
		}
		StsFile stsFile = StsFile.constructor(filename);
		StsAsciiFile asciiFile = new StsAsciiFile(stsFile);
		if (!asciiFile.openWrite())
		{
			new StsMessage(model.win3d, StsMessage.WARNING, "Failed to open file for writing: " + filename);
			return null;
		}
		// Break header into lines
		String exportHeader = getHeader() + getLineHeader(line);
		String textLine;
		while (exportHeader.indexOf("\n") != -1)
		{
			int index = exportHeader.indexOf("\n");
			textLine = StsStringUtils.padClipString(exportHeader.substring(0, index), 80);
			asciiFile.writeLine(textLine);
			exportHeader = exportHeader.substring(index + 1, exportHeader.length());
		}
		if (exportHeader.length() > 1)
			asciiFile.writeLine(exportHeader);
		return asciiFile;
	}

	public String getLineHeader(StsPreStackLine line)
    {
        return "vels for line: " + line.getName();
    }

    public Component getPanel()
	{
		return this;
	}
	public Component getPanel(boolean val) { return getPanel(); }

	public StsDialogFace getEditableCopy()
	{
		return (StsDialogFace)StsToolkit.copyObjectNonTransientFields(this);
	}

	public String getHeader()
	{
		if (textHeader == null)
			textHeader = getFileHeader();
		return textHeader;
	}

	protected void constructVelocityUnitsGroupBox()
	{
		unitsGroupBox = new StsGroupBox("Units and Scaling");
		velocityUnits = model.getProject().getVelocityUnits();
		unitsScale = 1.0f;
		StsComboBoxFieldBean unitsBean = new StsComboBoxFieldBean(this, "velocityUnitsString", "Velocity Units:", StsParameters.VEL_UNITS);
		unitsScaleBean = new StsFloatFieldBean(this, "unitsScale", false, "Units scale factor:");
		unitsGroupBox.addToRow(unitsBean);
		unitsGroupBox.addEndRow(unitsScaleBean);
		if(velocityUnits.equals(StsParameters.VEL_M_PER_MSEC) || velocityUnits.equals(StsParameters.VEL_M_PER_SEC)) 
		{
		    unitsBean.setValueObject(StsParameters.VEL_M_PER_SEC);
		}
		else
		{
		    unitsBean.setValueObject(StsParameters.VEL_FT_PER_SEC); 
		}
	}

	public void setExportName(String name)
	{
		this.exportName = name;
	}

	public String getExportName()
	{
		return exportName;
	}

	public void setNullValueString(String nullValueString)
	{
		for (int n = 0; n < 3; n++)
			if (nullValueStrings[n] == nullValueString)
				nullType = (byte)n;
	}

	public String getNullValueString()
	{
		return nullValueStrings[nullType];
	}

	public byte getNullType()
	{
		return nullType;
	}

	public void setVelocityTypeString(String type)
	{
		for (int n = 0; n < 4; n++)
			if (VEL_STRINGS[n] == type)
				velocityType = (byte)n;
	}

	public String getVelocityTypeString()
	{
		return VEL_STRINGS[velocityType];
	}

	public byte getVelocityType()
	{
		return velocityType;
	}

	public String getFileHeader()
	{
		if (exportHeader != null)
		{
			return exportHeader;
		}
		String hdr = "*** Hand Picked Velocity Profiles\n" +
		             "*** Exported from S2S Systems Software\n" +
		             "*** Date: " + new Date(System.currentTimeMillis()) + "\n";
		return hdr;
	}

	public float getRowNumMin()
	{
		return rowNumMin;
	}

	public float getRowNumMax()
	{
		return rowNumMax;
	}

	public float getColNumMin()
	{
		return colNumMin;
	}

	public float getColNumMax()
	{
		return colNumMax;
	}

	public float getZMin()
	{
		return zMin;
	}

	public float getZMax()
	{
		return zMax;
	}

	public float getRowNumInc()
	{
		return rowNumInc;
	}

	public float getColNumInc()
	{
		return colNumInc;
	}

	public float getZInc()
	{
		return zInc;
	}

	public void setRowNumMin(float value)
	{
		rowNumMin = value;
	}

	public void setRowNumMax(float value)
	{
		rowNumMax = value;
	}

	public void setRowNumInc(float value)
	{
		rowNumInc = value;
	}

	public void setColNumMin(float value)
	{
		colNumMin = value;
	}

	public void setColNumMax(float value)
	{
		colNumMax = value;
	}

	public void setColNumInc(float value)
	{
		colNumInc = value;
	}

	public void setZMin(float value)
	{
		zMin = value;
	}

	public void setZMax(float value)
	{
		zMax = value;
	}

	public void setZInc(float value)
	{
		zInc = value;
	}


	public void setExportDomainString(String string)
	{
		byte newExportDomainType = StsParameters.getZDomainFromString(string);
		if (newExportDomainType == exportDomainType)
			return;
		exportDomainType = newExportDomainType;
		zMin = volume.getZMin();
		zMax = volume.getZMax();
		zInc = volume.getZInc();
		if (exportDomainType == StsParameters.TD_DEPTH)
		{
			zMin = volume.getMinDepthAtTime(zMin);
			zMax = volume.getMaxDepthAtTime(zMax);
			int nSlices = volume.nSlices;
			double[] niceScale = StsMath.niceScale(zMin, zMax, nSlices, true);
			zMin = (float)niceScale[0];
			zMax = (float)niceScale[1];
			zInc = (float)niceScale[2];
            nSlices = Math.round((zMax - zMin)/zInc) + 1;
		}
	}

	public String getExportDomainString()
	{
		return StsParameters.getTDString(exportDomainType);
	}

	static public void main(String[] args)
	{
		StsPreStackVelocityModel3d velocityModel = new StsPreStackVelocityModel3d();
		StsModel model = new StsModel();
		StsProject project = new StsProject();
		model.setProject(project);
		StsPreStackHandVelExportPanel exportPanel = new StsPreStackHandVelExportPanel3d(model, velocityModel);
		StsOkCancelDialog exportDialog = new StsOkCancelDialog(null, exportPanel, "Test Export Panel", true);
	}

    protected abstract String getFileFooter();

    protected abstract String getOutputPathname(StsPreStackLine stsPreStackLine);

}
