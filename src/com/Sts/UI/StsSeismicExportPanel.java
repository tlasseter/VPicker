package com.Sts.UI;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;
import sun.io.*;

import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.nio.*;
import java.util.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */
public abstract class StsSeismicExportPanel extends StsJPanel implements StsDialogFace, StsSerializable
{
	protected StsModel model;
    protected StsProject project;
	protected StsSeismicBoundingBox volume = null;
    protected StsSeismicBoundingBox exportVolume;
    protected StsCroppedBoundingBox croppedBoundingBox;

	protected String panelTitle = defaultPanelTitle;
	protected boolean process = true;
    /** input volume geometry box */
    protected StsGroupBox volumeGroupBox;
    /** export volume geometry box */
    protected StsGroupBox exportGroupBox;
    /** stemname of exported file:  exportName.sgy */
	protected String exportName = null;
    /** user selected output sample format: IBM_FLOAT, IEEE, SCALED8, BYTE */
	protected byte outputFormat = IBM_FLOAT;
    /** identifies null value of input volume; these sampe values are converted to a corresponding output nul */
    protected byte nullType = NULL_ZERO;
    /** multiplier for samples on output */
	protected float sampleScale = 1.0f;
    /** flag indicated domain conversion from input to output: 0 = NONE, 01 = DEPTH->TIME, 10 = TIME->DEPTH, 03 = DEPTH->SEISMIC_DEPTH, 30 = SEISMIC_DEPTH->DEPTH */
    protected byte conversionFlag = StsParameters.TD_CONV_NONE;
    protected String exportHeader = null;
	//protected String filename;
    /** cdp x and y scalar in traceHeaders */
	protected float unitsScale = 1.0f;
    /** number of exported rows */
	protected int nRows;
	/** number of exported cols */
	protected int nCols;
	/** number of exported samples in trace */
	protected int nCroppedSlices;
	/** slice range for exported volume */
	protected int sliceStart, sliceEnd, sliceInc;
	/** number of rows in full volume */
	protected int nVolumeRows;
	/** number of cols in full volume */
	protected int nVolumeCols;
	/** number of samples in full trace */
	protected int nVolumeSlices;
	/** sample size of full volume */
	protected int sampleSize;
    /** float value output if sample is null */
	protected float nullFloat;
    /** scratch array for trace bytes */
    protected byte[] outputTraceBytes;
    /** scratch array for trace floats */
    protected float[] outputTraceFloats;
    /** domain (time, seismic depth, or depth) of input volume to export */
    protected String volumeDomainString = StsParameters.TD_NONE_STRING;
    /** domain (time, seismic depth, or depth) of exported output volume */
    protected String exportDomainString = StsParameters.TD_NONE_STRING;
    /** format of samples to be exported */
	protected String inputFormatString = inputFormatByte;
    /** vertical units of input volume to export */
    protected String volumeVerticalUnitsString;
    /** vertical units of output volume to export */
    protected String exportVerticalUnitsString;
    /** velocity volume export is handled differently */
	private byte volumeType;
	/** velocity type to export: rms, interval, instantaneous */
	protected byte velocityType = StsParameters.SAMPLE_TYPE_VEL_RMS;
	/** output units for velocity */
	protected String velocityUnits = StsParameters.VEL_FT_PER_MSEC;
    /** output vertical units for time or depth */
    protected String verticalUnits = StsParameters.D_NONE_STRING;
    /** velocity model to use if needed */
    protected StsSeismicVelocityModel velocityModel;

    protected float depthDatum = 0.0f;
	Frame frame;
    // these members are associated with the unitsGroupBox which is added to the seismicExportPanel
	protected StsGroupBox unitsGroupBox = null;
	protected StsFloatFieldBean unitsScaleBean;
	protected StsFloatFieldBean totalScaleBean;

	protected StsFile out = null;
	protected OutputStream os = null;
	protected BufferedOutputStream bos = null;
	protected ByteArrayOutputStream baos = null;
	protected DataOutputStream dos = null;

	protected byte[] textHeader = new byte[textHeaderSize];
	protected byte[] binaryHeader = new byte[binaryHeaderSize];
	protected byte[] traceHeader = new byte[traceHeaderSize];
	protected String stringTextHeader = null;

	JLabel fileViewLbl;

	protected StsGroupBox exportParametersGroupBox = new StsGroupBox("File Export Parameters");
	protected StsStringFieldBean exportNameBean;
	protected StsComboBoxFieldBean nullValueBean;
	// output format selection; if scaled float selected, sampleScale is enabled; text explanation displayed for each selection
	protected StsGroupBox formatGroupBox = new StsGroupBox("Select Output Format");
	protected StsStringFieldBean inputFormatBean;
	protected StsComboBoxFieldBean outputFormatBean;

    protected StsStringFieldBean volumeVerticalUnitsBean;
    protected StsComboBoxFieldBean exportVerticalUnitsBean;

    protected float traceHeaderScale = 1.0f;
	protected boolean isInputFloats = true;
	protected boolean isOutputFloats = true;

	protected boolean canceled = false;
	// null value selection
//	StsButtonListFieldBean nullValuesButtonListBean = new StsButtonListFieldBean(this, "nullValueString", "Set null values to:", nullValueStrings, true);

	// view text header
	StsButton editTextHeader = new StsButton("View/Edit Text Header", "View the file text header", this, "editTextHeader");
	static final String defaultPanelTitle = "SEGY PostStack3d Export";
	static final String defaultExportTypeString = ".segy";

    /** data is scaled bytes; output is IBM floats */
    static public final byte IBM_FLOAT = 0;
    /** data is scaled bytes; output is bytes */
	static public final byte BYTE = 1;
    /** data is scaled bytes; output is scaled bytes */
	static public final byte SCALED8 = 2;
    /** data is scaled bytes; output is IEEE_FLOAT */
	static public final byte IEEE_FLOAT = 3;

	static public final String BYTE_STRING = "Signed byte output format";
	static public final String SCALED8_STRING = "Scaled byte output format";
	static public final String IEEE_FLOAT_STRING = "IEEE Float output format";
	static public final String IBM_FLOAT_STRING = "IBM Float output format";
	static public final String[] outputFormatStrings = new String[] { IBM_FLOAT_STRING, BYTE_STRING, SCALED8_STRING, IEEE_FLOAT_STRING };
	static public final byte[] outputFormatBytes = new byte[] { IBM_FLOAT, BYTE, SCALED8, IEEE_FLOAT };
	static public final String inputFormatByte = "Scaled byte data format.";
	static public final String inputFormatFloat = "Float data format.";

	static public final byte NULL_NONE = 0;
	static public final byte NULL_MIN = 1;
	static public final byte NULL_MAX = 2;
	static public final byte NULL_ZERO = 3;

	static public final String NULL_NONE_STRING = "None";
	static public final String NULL_MIN_STRING = "Min";
	static public final String NULL_MAX_STRING = "Max";
	static public final String NULL_ZERO_STRING = "Zero";
	static public final String[] nullValueStrings = new String[] { NULL_NONE_STRING, NULL_MIN_STRING, NULL_MAX_STRING, NULL_ZERO_STRING };

	static public final String FLOAT_FORMAT_IEEE = "IEEE Float";
	static public final String FLOAT_FORMAT_IBM = "IBM Float";
	static public final String[] floatFormatStrings = new String[] { FLOAT_FORMAT_IEEE, FLOAT_FORMAT_IBM };

	static public final String[] VEL_STRINGS = StsParameters.VEL_STRINGS;

	static public final int textHeaderSize = StsSEGYFormat.defaultTextHeaderSize;
	static public final int binaryHeaderSize = StsSEGYFormat.defaultBinaryHeaderSize;
	static public final int traceHeaderSize = StsSEGYFormat.defaultTraceHeaderSize;

	static public final byte nullByte = StsParameters.nullByte;

	static final String group = "velModelExport";
	static final String format = "txt";

	abstract protected void doExportSeismic(StsProgressPanel progressPanel);
    abstract protected StsGroupBox getInputGeometryGroupBox();
    abstract protected StsGroupBox getOutputGeometryGroupBox();

	public StsSeismicExportPanel()
	{
	}

	public StsSeismicExportPanel(StsModel model, StsSeismicBoundingBox volume, String title, boolean isAVelocityExport)
	{
		this.volumeType = volume.volumeType;
		this.model = model;
        this.project = model.getProject();
		this.volume = volume;
		exportName = StsStringUtils.cleanString(volume.getName() + "export");
		this.panelTitle = title;
		volume.setAngle();
        exportVolume = new StsSeismicBoundingBox(false);
        StsToolkit.copyAllObjectFields(volume, exportVolume, true);
		frame = model.win3d;
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
		fileViewLbl = new JLabel(panelTitle);

		exportNameBean = new StsStringFieldBean(this, "exportName", "Export name");
		exportNameBean.setToolTipText("Exports to file directory/name.segy");

		inputFormatBean = new StsStringFieldBean(this, "inputFormatString", false, "Current data format:");
		outputFormatBean = new StsComboBoxFieldBean(this, "outputFormatString", "Output Format", outputFormatStrings);
		nullValueBean = new StsComboBoxFieldBean(this, "nullValueString", "Null value:", this.nullValueStrings);

		velocityType = StsParameters.SAMPLE_TYPE_VEL_RMS;
	}

	protected void assembleGroupBoxes()
	{
		// constructVelocityUnitsGroupBox();

		fileViewLbl.setFont(new java.awt.Font("Serif", 1, 14));
        exportParametersGroupBox.gbc.fill = HORIZONTAL;
		exportParametersGroupBox.addToRow(exportNameBean);

		formatGroupBox.gbc.gridwidth = 2;
        formatGroupBox.gbc.fill = HORIZONTAL;
		formatGroupBox.add(inputFormatBean);
		formatGroupBox.add(outputFormatBean);
		formatGroupBox.add(nullValueBean);
     /*
		if (isVelocity())
		{
			StsParameterFileFieldBean templateBean = new StsParameterFileFieldBean(this, StsSeismicExportPanel.class, StsSeismicExportPanel.class, group, format);
			templateBean.setComboBoxToolTipText("Select existing velocity model export template");
			// templateBean.setNewButtonToolTipText("Create a new velocity model export template and save it");
			StsComboBoxFieldBean velocityTypeBean = new StsComboBoxFieldBean(this, "velocityTypeString", "Velocity type", StsParameters.VEL_STRINGS);
			exportParametersGroupBox.addToRow(velocityTypeBean);
			// exportParametersGroupBox.addToRow(templateBean, 2, 1.0);
		}
    */
	}

	private void constructVelocityUnitsGroupBox()
	{
		unitsGroupBox = new StsGroupBox("Output Velocity Units and Scaling");
		velocityUnits = model.getProject().getVelocityUnits();
        verticalUnits = model.getProject().getVerticalUnitsString();
		unitsScale = 1.0f;
		StsComboBoxFieldBean unitsBean = new StsComboBoxFieldBean(this, "velocityUnitsString", "Velocity Units:", StsParameters.VEL_UNITS);

		unitsScaleBean = new StsFloatFieldBean(this, "unitsScale", false, "Units scale factor:");
		unitsGroupBox.addToRow(unitsBean);
		unitsGroupBox.addEndRow(unitsScaleBean);
	}

	void assemblePanel()
	{
		gbc.fill = GridBagConstraints.HORIZONTAL;
		add(fileViewLbl);
		add(exportParametersGroupBox);
		//add(subVolumeGroupBox);
		volumeGroupBox = getInputGeometryGroupBox();
        add(volumeGroupBox);
        exportGroupBox = getOutputGeometryGroupBox();
        add(exportGroupBox);
		add(formatGroupBox);
		if (unitsGroupBox != null)
			add(unitsGroupBox);
		gbc.fill = GridBagConstraints.NONE;
		add(editTextHeader);
		gbc.fill = GridBagConstraints.HORIZONTAL;
		//add(progressPanel);
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

	public float getUnitsScale()
	{
		return unitsScale;
	}

	public void setUnitsScale(float scale)
	{
		this.unitsScale = scale;
	}

	public String getVelocityUnits()
	{
		return velocityUnits;
	}

	public void setVelocityTypeString(String type)
	{
		for(int n = 0; n < 4; n++)
		{
			if (VEL_STRINGS[n] == type)
				velocityType = (byte)n;
		}
		volume.setVelocityTypeString(type);
	}

	public String getVelocityTypeString()
	{
		return VEL_STRINGS[velocityType];
	}

	public byte getVelocityType()
	{
		return velocityType;
	}

    public void setExportName(String name)
    {
        this.exportName = name;
    }

    public String getExportName()
    {
        return exportName;
    }

	public void setOutputFormatString(String outputString)
	{
		outputFormat = IBM_FLOAT;
		for (int n = 0; n < 4; n++)
		{
			if (outputFormatStrings[n] == outputString)
			{
				outputFormat = outputFormatBytes[n];
				break;
			}
		}
//		typeDesc.setText(getFormatString());
		if (outputFormat == BYTE)
			sampleScale = 1.0f;
	}

	public String getOutputFormatString()
	{
		return outputFormatStrings[outputFormat];
	}


    public String getInputFormatString()
    {
        return inputFormatString;
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

	public void setSampleScale(float scale)
	{
		this.sampleScale = scale;
	}

	public float getSampleScale()
	{
		return sampleScale;
	}

	public void editTextHeader()
	{
		constructTextHeader();

		StsTextHeaderEditDialog volumeDialog = new StsTextHeaderEditDialog(frame, stringTextHeader, false);
		volumeDialog.setVisible(true);
		if (volumeDialog.getMode() == volumeDialog.OK)
		{
			stringTextHeader = volumeDialog.getHeaderText();
			textHeader = stringTextHeader.getBytes();
			textHeader = (byte[])StsMath.trimArray(textHeader, textHeaderSize);
		}
	}

	public void dialogSelectionType(int type)
	{
		if (type == StsDialogFace.PROCESS || type == StsDialogFace.OK)
		{
			run();
		}
		if (type == StsDialogFace.CANCEL)
		{
			canceled = true;
		}
	}
	public Component getPanel(boolean val) { return getPanel(); }
	public Component getPanel()
	{
		return this;
	}

	public StsDialogFace getEditableCopy()
	{
		return (StsDialogFace)StsToolkit.copyObjectNonTransientFields(this);
	}

	public boolean run()
	{
		if (!checkExportSeismic())
		{
			return false;
		}
		final StsDialogFace face = this;
		Runnable runExport = new Runnable()
		{
			public void run()
			{
				StsProgressCancelDialog dialog = null;
				try
				{
					dialog = new StsProgressCancelDialog(model.win3d, face, "Export Seismic Volume", "Export Progress", false);
					StsProgressPanel progressPanel = dialog.getProgressPanel();
					dialog.setVisible(true);
					doExportSeismic(progressPanel);
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

	abstract protected boolean checkExportSeismic();

	abstract protected boolean setupExport();

	protected void setupFile(StsSEGYFormat segyFmt, String fileName) throws IOException
	{
		out = StsFile.constructor(fileName);
		dos = out.getDataOutputStream();
	}

    protected void computeConversionFlag()
    {
        byte volumeZDomain = volume.getZDomain();
        byte exportZDomain = exportVolume.getZDomain();
        conversionFlag = StsParameters.computeDomainConversionFlag(volumeZDomain, exportZDomain);
        if(conversionFlag == StsParameters.TD_CONV_NONE) return;

        velocityModel = model.getProject().velocityModel;
        if(velocityModel == null)
            conversionFlag = StsParameters.TD_CONV_NONE;
    }

	protected void convertInputBytesToOutputBytes(byte[] inputBytes, byte[] outputBytes, int offset, int nSlices)
	{
		if (isVelocity())
		{
			convertInputBytesToOutputBytesVelocity(inputBytes, outputBytes, offset, nSlices);
		}
		else
		{
			convertInputBytesToOutputBytesVolume(inputBytes, outputBytes, offset, nSlices);
		}
	}

    protected boolean isVelocity() { return StsParameters.volumeTypeIsVelocity(volumeType); }

    /** input bytes are "C" unsigned bytes (0-254 with 255 as the null), but in Java which is signed, they will range
     *  from 0-127 then -128 to -2 with -1 as the null (0->0, 127->127, -128->128, -2->254, -1->255).
     *  On output, we want to convert these unsigned bytes to "C" signed bytes.
     *  Since Java doesn't support unsigned bytes, we are converting from signed equivalents.
     *  unsigned    JavaSignedEquivalent  signedResult
     *    0              0                -128
     *    127            127              -1
     *    128            -128             0
     *    254            -2               126
     *    255            -1               127
     * @param inputBytes
     * @param outputBytes
     * @param offset
     * @param nSlices
     */
    private void convertInputBytesToOutputBytesVolume(byte[] inputBytes, byte[] outputBytes, int offset, int nSlices)
	{
		float z = exportVolume.zMin;
		for (int n = 0; n < nSlices; n++, z += exportVolume.zInc)
		{
			int sliceNum = Math.round((z - exportVolume.zMin) / exportVolume.zInc);
			outputBytes[n] = StsMath.unsignedByteToSignedByte(inputBytes[offset + sliceNum]);
		}
	}

	private void convertInputBytesToOutputBytesVelocity(byte[] inputBytes, byte[] outputBytes, int offset, int nSlices)
	{
		if(exportDomainString == StsParameters.TD_TIME_STRING)
			convertInputBytesToOutputBytesVolume(inputBytes, outputBytes, offset, nSlices);
		else
		{
			double t1 = volume.zMin;
			double v1 = volume.getScaledValue(inputBytes[offset])/2;
			double z1 = depthDatum;
			double z = exportVolume.zMin;
			double dt = exportVolume.zInc;
			int n = 0;
			byte vInt = 0;
			for(int i = 1; i < nVolumeSlices; i++)
			{
				double v0 = v1;
				double t0 = t1;
				v1 = volume.getScaledValue(inputBytes[offset + i])/2;
				t1 += exportVolume.zInc;
				z1 = v1*t1 + depthDatum;
				double dv = v1 - v0;
				if(velocityType == StsParameters.SAMPLE_TYPE_VEL_INTERVAL)
				{
					float vIntF = (float)Math.sqrt((v1 * v1 * t1 - v0 * v0 * t0) / exportVolume.zInc);
					vInt = volume.getByteValue(vIntF);
				}
				while(z <= z1)
				{
					if(velocityType == StsParameters.SAMPLE_TYPE_VEL_RMS)
					{
						double f = (z - v0 * t0 - depthDatum) / (t0 * dv + v0 * dt);
						outputBytes[n++] = volume.getByteValue(v0 + f * dv);
					}
					else
						outputBytes[n++] = vInt;

					if(n == nVolumeSlices) return;
					z += exportVolume.zInc;
				}
			}
		}
	}

	protected void convertOutputBytesToOutputFloats(byte[] outputBytes, float[] outputFloats)
	{
		for (int n = 0; n < nVolumeSlices; n++)
		{
			byte outputByte = outputBytes[n];
			if (outputByte == nullByte)
				outputFloats[n] = nullFloat;
			else
				outputFloats[n] = volume.getFloatFromSignedByte(outputByte) * sampleScale;
		}
	}

	protected boolean writeOutputFloats(float[] outputFloats, int nSlices)
	{
		if (outputFormat == IBM_FLOAT)
		{
			for (int n = 0; n < nSlices; n++)
				if (!writeIbmFloat(outputFloats[n] * unitsScale))
					return false;
		}
		else
		{
			try
			{
				for (int n = 0; n < nSlices; n++)
					dos.writeFloat(outputFloats[n] * unitsScale);
			}
			catch (Exception e)
			{
				StsException.outputException("StsSeismicExportPanel.writeOutputFloats() failed.", e, StsException.WARNING);
				return false;
			}
		}
		return true;
	}

	protected boolean writeOutputBytes(byte[] outputBytes)
	{
		try
		{
            dos.write(outputBytes);
		}
		catch (Exception e)
		{
			StsException.outputException("StsSeismicExportPanel.writeOutputFloats() failed.", e, StsException.WARNING);
			return false;
		}
		return true;
	}

	protected void convertInputFloatsToOutputFloats(float[] inputFloats, float[] outputFloats, int offset, int nSlices, int row, int col)
	{
		if (isVelocity())
		{
            if(conversionFlag == StsParameters.TD_CONV_NONE)
			    convertInputFloatsToOutputFloatsVolume(inputFloats, outputFloats, offset, nSlices);
            else if(conversionFlag == StsParameters.TD_CONV_TIME_DEPTH)
            {
                float x = volume.getXCoor(col);
                float y = volume.getYCoor(row);
                convertInputFloatsToOutputFloatsVelocityTimeToDepth(inputFloats, outputFloats, offset, nSlices, x, y);
            }
            else if(conversionFlag == StsParameters.TD_CONV_DEPTH_TIME)
            {
                float x = volume.getXCoor(col);
                float y = volume.getYCoor(row);
                convertInputFloatsToOutputFloatsVelocityDepthToTime(inputFloats, outputFloats, offset, nSlices, x, y);
            }
        }
		else
		{
            if(conversionFlag == StsParameters.TD_CONV_NONE)
			    convertInputFloatsToOutputFloatsVolume(inputFloats, outputFloats, offset, nSlices);
            else if(conversionFlag == StsParameters.TD_CONV_TIME_DEPTH)
            {
                float x = volume.getXCoor(col);
                float y = volume.getYCoor(row);
                convertInputFloatsToOutputFloatsVolumeTimeToDepth(inputFloats, outputFloats, offset, nSlices, x, y);
            }
            else if(conversionFlag == StsParameters.TD_CONV_DEPTH_TIME)
            {
                float x = volume.getXCoor(col);
                float y = volume.getYCoor(row);
                convertInputFloatsToOutputFloatsVolumeDepthToTime(inputFloats, outputFloats, offset, nSlices, x, y);
            }
        }
	}

	private void convertInputFloatsToOutputFloatsVolume(float[] inputFloats, float[] outputFloats, int offset, int nSlices)
	{
		float z = croppedBoundingBox.zMin;
        float zInc = croppedBoundingBox.zInc;
		for (int n = 0; n < nSlices; n++, z += zInc)
		{
			int sliceNum = Math.round((z - volume.zMin) / volume.zInc);
			outputFloats[n] = inputFloats[offset + sliceNum];
		}
	}

	private void convertInputFloatsToOutputFloatsVolumeTimeToDepth(float[] inputFloats, float[] outputFloats, int offset, int nSlices, float x, float y)
	{
		float z = croppedBoundingBox.zMin;
        float test = volume.zMin;
		for (int n = 0; n < nSlices; n++, z += croppedBoundingBox.zInc)
		{
            float t = (float)velocityModel.getT(x, y, z, test);
            test = t;
            int sliceNum = volume.getNearestBoundedSliceCoor(t);
			outputFloats[n] = inputFloats[offset + sliceNum];
		}
	}

	private void convertInputFloatsToOutputFloatsVolumeDepthToTime(float[] inputFloats, float[] outputFloats, int offset, int nSlices, float x, float y)
	{
        float t = croppedBoundingBox.zMin;
        float tInc = croppedBoundingBox.zInc;
		for (int n = 0; n < nSlices; n++, t += tInc)
		{
            float z = (float)velocityModel.getZ(x, y, t);
            int sliceNum = volume.getNearestBoundedSliceCoor(z);
			outputFloats[n] = inputFloats[offset + sliceNum];
		}
	}

    // TODO implement TimeToTime and DepthToDepth going between AVG_VEL and INT_VEL
	private void convertInputFloatsToOutputFloatsVelocity(float[] inputFloats, float[] outputFloats, int offset, int nSlices)
	{
        if(volume.volumeType == exportVolume.volumeType)
            convertInputFloatsToOutputFloatsVolume(inputFloats, outputFloats, offset, nSlices);
	}

	private void convertInputFloatsToOutputFloatsVelocityTimeToDepth(float[] inputFloats, float[] outputFloats, int offset, int nSlices, float x, float y)
	{
        if(volume.volumeType == exportVolume.volumeType)
            this.convertInputFloatsToOutputFloatsVolumeTimeToDepth(inputFloats, outputFloats, offset, nSlices, x, y);
        else if(volume.volumeType == StsParameters.SAMPLE_TYPE_VEL_AVG)  // exportVolume.volumeType == SAMPLE_TYPE_VEL_INTERVAL
        {
            float tInc = volume.zInc;
            float zMin = croppedBoundingBox.zMin;
            float zInc = croppedBoundingBox.zInc;
            float z = zMin;
            float test = volume.zMin;
            for (int n = 0; n < nSlices; n++, z += zInc)
            {
                float t0 = (float)velocityModel.getT(x, y, z, test);
                test = t0;
                int sliceNum = volume.getBelowBoundedSliceCoor(t0);
                float v0 = inputFloats[offset + sliceNum];
                float v1 = inputFloats[offset + sliceNum + 1];
                float t1 = t0 + volume.zInc;
                outputFloats[n] = (float)Math.sqrt((v1 * v1 * t1 - v0 * v0 * t0) /tInc);
            }
        }
        else if(volume.volumeType == StsParameters.SAMPLE_TYPE_VEL_INTERVAL)  // exportVolume.volumeType ==  SAMPLE_TYPE_VEL_AVG
        {
            float tInc = volume.zInc;
            float zMin = croppedBoundingBox.zMin;
            float zInc = croppedBoundingBox.zInc;
            float z = zMin;
            float test = volume.zMin;
            float depthDatum = velocityModel.depthDatum;
            float timeDatum = velocityModel.timeDatum;
            for (int n = 0; n < nSlices; n++, z += zInc)
            {
                float t0 = (float)velocityModel.getT(x, y, z, test);
                outputFloats[n] = (z - depthDatum)/(t0 - timeDatum);
            }
        }
	}

	private void convertInputFloatsToOutputFloatsVelocityDepthToTime(float[] inputFloats, float[] outputFloats, int offset, int nSlices, float x, float y)
	{
	}

	protected void convertOutputFloatsToOutputBytes(float[] outputFloats, byte[] outputBytes, int nSlices)
	{
		for (int n = 0; n < nSlices; n++)
		{
			float outputFloat = outputFloats[n];
			if (outputFloat == nullFloat)
				outputBytes[n] = nullByte;
			else
				outputBytes[n] = volume.getByteValue(outputFloat);
		}
	}

	final boolean writeIbmFloat(float value)
	{
		int ibmBits = StsParameters.nullInteger;
		try
		{
			ibmBits = StsMath.convertIeeeFloatToIbmBits(value, false);
		}
		catch (Exception e)
		{
			new StsMessage(frame, StsMessage.WARNING, "Failed to convert float " + value + " to IBM float format.");
			return false;
		}
		try
		{
			dos.writeInt(ibmBits);
		}
		catch (Exception e)
		{
			new StsMessage(frame, StsMessage.WARNING, "Failed to write IBM float " + value + " to output file.");
			return false;
		}
		return true;
	}

	public void constructTextHeader()
	{
		StringBuffer stringBuffer = new StringBuffer(textHeaderSize);
		for (int n = 0; n < textHeaderSize; n++)
			stringBuffer.append(' ');

		int length = 80;
		int nLine = 0;

		putString(stringBuffer, nLine++, length, "-Created from S2S Systems Software");
        putString(stringBuffer, nLine++, length, "-Generated from " + StsToolkit.getSimpleClassname(volume) + " " + volume.getName());
        putString(stringBuffer, nLine++, length, "TraceNumberHeaderByte 1");
		putString(stringBuffer, nLine++, length, "InLineHeader Byte 9");
		putString(stringBuffer, nLine++, length, "CrosslineHeaderByte 13 ");
		putString(stringBuffer, nLine++, length, "CoordinateScalarByte 71");
		putString(stringBuffer, nLine++, length, "XCoordinate Byte 73");
		putString(stringBuffer, nLine++, length, "YCoordinateByte 77");
		putString(stringBuffer, nLine++, length, "NumberOfSamples Byte 115");
		putString(stringBuffer, nLine++, length, "SampleIntervalByte 117");
		putString(stringBuffer, nLine++, length, "DataScalarByte 169");
        putString(stringBuffer, nLine++, length, "Domain " + exportVolume.zDomain);
        putString(stringBuffer, nLine++, length, "SampleType " + StsParameters.getSampleTypeString(volume.volumeType));
        putString(stringBuffer, nLine++, length, "SampleUnits " + volume.getUnitsString());
		double[] origin = volume.getAbsoluteXY(croppedBoundingBox.xMin, croppedBoundingBox.yMin);
		putString(stringBuffer, nLine++, length, "XOrigin " + origin[0]);
        putString(stringBuffer, nLine++, length, "YOrigin " + origin[1]);
        putString(stringBuffer, nLine++, length, "Angle " + volume.angle);
        putString(stringBuffer, nLine++, length, "XSize " + croppedBoundingBox.getXSize());
        putString(stringBuffer, nLine++, length, "YSize " + croppedBoundingBox.getYSize());
        putString(stringBuffer, nLine++, length, "ZMin " + croppedBoundingBox.zMin);
        putString(stringBuffer, nLine++, length, "ZMax " + croppedBoundingBox.zMax);
        putString(stringBuffer, nLine++, length, "ZInc " + croppedBoundingBox.zInc);
		putString(stringBuffer, nLine++, length, getOutputFormatString() + " " + outputFormat);

        for(; nLine < 40; nLine++)
            putString(stringBuffer, nLine, length, "");

		stringTextHeader = stringBuffer.toString();
		textHeader = stringTextHeader.getBytes();
		textHeader = (byte[])StsMath.trimArray(textHeader, textHeaderSize);
	}

	static private void putString(StringBuffer stringBuffer, int nLine, int length, String string)
	{
        if(nLine < 9)
            string = "C" + (nLine+1) + "  " + string;
        else
            string = "C" + (nLine+1) + " " + string;
		char[] chars = string.toCharArray();
		int stringLength = Math.min(chars.length, length);
		stringBuffer.insert(nLine*length, chars, 0, stringLength);
	}

	static private void putString(ByteBuffer byteBuffer, String string, int offset, int length)
	{
		byte[] bytes = string.getBytes();
		length = Math.min(bytes.length, length);
		byteBuffer.put(bytes, offset, length);
	}

	protected void writeTextHdr(StsSEGYFormat fmt)
	{
		CharToByteCp500 converter = new CharToByteCp500();
		try
		{
			char[] chars = stringTextHeader.toCharArray();
			byte[] bytes = converter.convertAll(chars);
			bytes = (byte[])StsMath.trimArray(bytes, textHeaderSize);
			dos.write(bytes);
		}
		catch (Exception e)
		{
			return;
		}
	}

    protected void writeBinaryHeader()
    {
		try
		{
            dos.write(binaryHeader);
		}
		catch (Exception e)
		{
			return;
		}
    }

    public void constructTraceHdr(StsSEGYFormat fmt, double il, double xl, double x, double y, double ns, double si,
											double type, double tn, double cs, double scale)
	{
		try
		{
			ByteBuffer wrappedBytes = ByteBuffer.wrap(traceHeader);
			wrappedBytes.putInt(0, (int)tn);
			wrappedBytes.putInt(8, (int)il);
			wrappedBytes.putInt(12, (int)xl);
			wrappedBytes.putShort(28, (short)type);
			wrappedBytes.putShort(70, (short)cs);
			wrappedBytes.putInt(72, (int)x);
			wrappedBytes.putInt(76, (int)y);
			wrappedBytes.putShort(114, (short)ns);
			wrappedBytes.putShort(116, (short)si);
			wrappedBytes.putShort(168, (short)scale);
		}
		catch (Exception e)
		{
			new StsMessage(model.win3d, StsMessage.ERROR,
								"Failed writing header for inline " + il + " crossline " + xl + " -- " + e);
		}
	}

	protected void constructBinaryHdr()
	{
		try
		{
			ByteBuffer wrappedBytes = ByteBuffer.wrap(binaryHeader);
			wrappedBytes.putShort(16, (short)(croppedBoundingBox.zInc * 1000));
			wrappedBytes.putShort(20, (short) nCroppedSlices);
			if (outputFormat == BYTE || outputFormat == SCALED8)
				wrappedBytes.putShort(24, (short)8);
			else if (outputFormat == IEEE_FLOAT)
				wrappedBytes.putShort(24, (short)5);
			else if (outputFormat == IBM_FLOAT)
				wrappedBytes.putShort(24, (short)1);
		}
		catch (Exception e)
		{
			new StsMessage(model.win3d, StsMessage.ERROR, "Failed writing binary header -- " + e);
		}
	}

	static public void createDialog(StsModel model, StsSeismicBoundingBox volume, String title, boolean isAVelocityExport)
	{
		if (!volume.isArchived)
		{
			StsSeismicExportPanel exportPanel = null;
			if (volume instanceof StsPreStackVelocityModel2d)
			{
				exportPanel = new StsSeismicExportPanel2d(model, volume, title, isAVelocityExport);
			}
			else
			{
				exportPanel = new StsSeismicExportPanel3d(model, volume, title, isAVelocityExport);
			}
			new StsProcessDismissDialog(model.win3d, exportPanel, exportPanel, title, true);
		}
	}

	public String getFileHeader()
	{
		if (exportHeader != null)
		{
			return exportHeader;
		}
		String hdr = new String("Hand Picked Velocity Profiles\n" +
										"Exported from S2S Systems Software\n" +
										"Date: " + new Date(System.currentTimeMillis()) + "\n");
		return hdr;
	}

	public static void main(String[] args)
	{
		/*
			StringBuffer stringBuffer = new StringBuffer(160);
			for(int n = 0; n < 160; n++)
		 stringBuffer.append(' ');
			String string;
			string = "C1";
			putString(stringBuffer, 0, 80, string);
			string = "C2";
			putString(stringBuffer, 80, 80, string);
			String outputString = stringBuffer.toString();
			System.out.println(outputString);
			byte[] bytes = outputString.getBytes();
			bytes = (byte[])StsMath.trimArray(bytes, 160);
			System.out.println("bytes.length " + bytes.length);
		 */
		StsModel model = new StsModel();
		StsProject project = new StsProject();
		model.setProject(project);
		StsSeismicVolume volume = new StsSeismicVolume();
		volume.setName("testVolume");
		StsGroupBox extraGroupBox = new StsGroupBox("Extra");
		StsStringFieldBean extraBean = new StsStringFieldBean("Extra field.");
		extraGroupBox.add(extraBean);
		StsSeismicExportPanel.createDialog(model, volume, "Seismic Export Test", false);
	}
}
