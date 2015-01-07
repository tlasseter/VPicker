package com.Sts.Actions.Wizards.Seismic;

import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.lang.reflect.*;
import java.text.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2006</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

public class StsFileFormatPanel extends StsJPanel
{
//    private boolean isSampleFormatOverride = false;
    private String overrideSampleFormatString = StsSEGYFormat.NONE_STRING;
    private String headerSampleFormatString = StsSEGYFormat.IBMFLT_STRING;
    private boolean isLittleEndian = false;
    private int binaryHeaderSize = StsSEGYFormat.defaultBinaryHeaderSize;
    private int traceHeaderSize = StsSEGYFormat.defaultTraceHeaderSize;
    private String textHeaderFormatString = StsSEGYFormat.defaultTextHeaderFormatString;
    private int textHeaderSize = StsSEGYFormat.defaultTextHeaderSize;
    private float startZ = 0.0f;
    private boolean overrideHeader = false;
    private int overrideNSamples = 0;
    private float overrideSampleSpacing = 4.0f;

//    private StsBooleanFieldBean overrideSampleFormatBean = new StsBooleanFieldBean(this, "isSampleFormatOverride", "Override sample format", true);
    private StsStringFieldBean headerSampleFormatBean = new StsStringFieldBean(this, "headerSampleFormatString", false, "Header sample format:");
    private StsComboBoxFieldBean overrideSampleFormatComboBoxBean = new StsComboBoxFieldBean(this, "overrideSampleFormatString", "Override sample format:", StsSEGYFormat.sampleFormatStrings);
    private StsBooleanFieldBean endianCheckBoxBean = new StsBooleanFieldBean(this, "isLittleEndian", "Little Endian", true);
    private StsIntFieldBean binaryHeaderSizeBean = new StsIntFieldBean(this, "binaryHeaderSize", true, "Binary Header Size:");
    private StsIntFieldBean textHeaderSizeBean = new StsIntFieldBean(this, "textHeaderSize", true, "Text Header Size:");
    private StsComboBoxFieldBean textFmtComboBean = new StsComboBoxFieldBean(this, "textHeaderFormatString", "Text Header Format:", StsSEGYFormat.textHeaderFormatStrings);
    private StsIntFieldBean traceHeaderSizeBean = new StsIntFieldBean(this, "traceHeaderSize", 0, 960, "Trace Header Size:");
    private StsFloatFieldBean startZBean = new StsFloatFieldBean(this, "startZ", true, "Start Z:");
    private StsBooleanFieldBean overrideHdrBean = new StsBooleanFieldBean(this, "override", "Override:", true);
    private StsIntFieldBean overrideNSamplesBean = new StsIntFieldBean(this, "overrideNSamples", false, "Samples/Trace:");
    private StsFloatFieldBean overrideSampleSpacingBean = new StsFloatFieldBean(this, "overrideSampleSpacing", false, "Sample Spacing:");

	private StsButtonFieldBean viewBinHeaderButton = new StsButtonFieldBean("View Volume(s) Binary Header", "View the binary headers.", this, "viewBinaryHeaders");
	private StsButtonFieldBean viewTxtHeaderButton = new StsButtonFieldBean("View Volume(s) Text Header", "View the text headers.", this, "viewTextHeaders");

    {
        headerSampleFormatBean.getTextField().setColumns(10);
        textHeaderSizeBean.setColumns(4);
        textFmtComboBean.getComboBox().setMaximumSize(new Dimension(30, 10));
        binaryHeaderSizeBean.setColumns(4);
        traceHeaderSizeBean.setColumns(3);

        binaryHeaderSizeBean.setToolTipText("Specify the binary header size");
        textHeaderSizeBean.setToolTipText("Specify the text header size");
        textFmtComboBean.setToolTipText("Select the text header format");
        startZBean.setToolTipText("+ indicates above Sea Level");
        overrideHdrBean.setToolTipText("Disregard headers and allow user to supply samples/trace and sample spacing");
    }

    StsGroupBox sampleFormatBox = new StsGroupBox("Data Format");
    StsJPanel headerBox = new StsGroupBox("File and trace headers");
    StsJPanel headerSizePanel = new StsGroupBox();
    StsGroupBox traceGeometryBox = new StsGroupBox("Trace Geometry");
    StsSEGYFormat segyFormat = null;

	private StsSeismicWizard wizard = null;
	private boolean editable = true;

    protected Method get = null;
	protected Method set = null;

    JTextField messageTxt = new JTextField();

    private DecimalFormat fmt = new DecimalFormat("#####");
    private static final int LINE_LENGTH = 80;

    public StsFileFormatPanel(StsSeismicWizard seismicWizard, boolean editable)
	{
		this.editable = editable;
		this.wizard = seismicWizard;
		constructPanel();			
	}

    /** We have a single segyFormat which is user defined or the default.  A few parameters (sampleFormatOverride and sampleFormatvalue)
     *  are extracted from it and isVisible here.  If these or other parameters are changed, they are applied to the selectedVolumes
     *  (specifically the segyData members of the selectedVolumes).  We this wizardStep is finished, we set the two segyFormat parameters
     *  back in segyFormat where they will be saved (if changed) in a user-defined segyFormat.
     */
    public void initialize(StsSegyData[] segyDatasets)
    {
        segyFormat = wizard.getSegyFormat();

//        isSampleFormatOverride = segyFormat.isSampleFormatOverride;
 //       overrideSampleFormatBean.setValue(isSampleFormatOverride);

//        sampleFormatComboBoxBean.setEditable(isSampleFormatOverride);

        int sampleFormat = segyDatasets[0].getHeaderSampleFormat();
        for(int n = 1; n < segyDatasets.length; n++)
        {
            int otherSampleFormat = segyDatasets[n].getHeaderSampleFormat();
            if(otherSampleFormat != sampleFormat)
            {
                sampleFormat = StsSEGYFormat.NONE;
                break;
            }
        }
        headerSampleFormatBean.setValue(StsSEGYFormat.sampleFormatStrings[sampleFormat]);

        sampleFormat = segyDatasets[0].getOverrideSampleFormat();
        for(int n = 1; n < segyDatasets.length; n++)
        {
            int otherSampleFormat = segyDatasets[n].getOverrideSampleFormat();
            if(otherSampleFormat != sampleFormat)
            {
                sampleFormat = StsSEGYFormat.NONE;
                break;
            }
        }
        overrideSampleFormatComboBoxBean.doSetValueObject(StsSEGYFormat.sampleFormatStrings[sampleFormat]);

        endianCheckBoxBean.setValue(segyDatasets[0].isLittleEndian);
		binaryHeaderSizeBean.setValue(segyDatasets[0].binaryHeaderSize);
		textHeaderSizeBean.setValue(segyDatasets[0].textHeaderSize);
        traceHeaderSizeBean.setValue(segyDatasets[0].traceHeaderSize);

        if(!overrideHeader)
        {
            overrideSampleSpacingBean.setValue(segyDatasets[0].getSampleSpacing());
            overrideNSamplesBean.setValue(segyDatasets[0].getNSamples());
        }
    }

	private void constructPanel()
	{
        gbc.fill = gbc.HORIZONTAL;

        sampleFormatBox.gbc.fill = gbc.HORIZONTAL;
        sampleFormatBox.addToRow(headerSampleFormatBean, 2, 1.0);
        sampleFormatBox.addToRow(overrideSampleFormatComboBoxBean, 2, 1.0);
        sampleFormatBox.addEndRow(endianCheckBoxBean);
        add(sampleFormatBox);

        headerBox.gbc.fill = gbc.HORIZONTAL;

        StsJPanel headerSizePanel = StsJPanel.noInsets();
        headerSizePanel.gbc.fill = gbc.HORIZONTAL;
        headerSizePanel.addToRow(textHeaderSizeBean, 2, 1.0);
		headerSizePanel.addToRow(textFmtComboBean, 2, 1.0);
        headerSizePanel.addToRow(binaryHeaderSizeBean);
		headerSizePanel.addEndRow(traceHeaderSizeBean);
        headerBox.add(headerSizePanel);

        StsJPanel viewHeadersPanel = StsJPanel.noInsets();
        viewHeadersPanel.gbc.fill = gbc.HORIZONTAL;
        viewHeadersPanel.addToRow(viewBinHeaderButton);
		viewHeadersPanel.addToRow(viewTxtHeaderButton);
        headerBox.add(viewHeadersPanel);
        
        add(headerBox);

        traceGeometryBox.gbc.gridx = 1;
        traceGeometryBox.gbc.fill = gbc.HORIZONTAL;
        traceGeometryBox.addToRow(startZBean);
        traceGeometryBox.addToRow(overrideHdrBean);
		traceGeometryBox.addToRow(overrideNSamplesBean);
		traceGeometryBox.addEndRow(overrideSampleSpacingBean);
		add(traceGeometryBox);
	}

    public void updatePanel()
    {
        wizard.setSkipReanalyzeTraces(true);
        StsSeismicBoundingBox[] selectedVolumes = wizard.getSelectedVolumes();
        StsSegyData[] selectedSegyDatasets = StsSegyData.getSegyDatasets(selectedVolumes);
        initialize(selectedSegyDatasets);
        wizard.setSkipReanalyzeTraces(false);
//        wizard.analyzeHeaders();
    }
/*
    public void setIsSampleFormatOverride(boolean b)
    {
        if(isSampleFormatOverride == b) return;
        isSampleFormatOverride = b;
        sampleFormatComboBoxBean.setEditable(b);
        segyFormat.setIsSampleFormatOverride(b);
        StsSeismicBoundingBox[] vols = wizard.getSelectedVolumes();
    	for(int i=0; i<vols.length; i++)
    		vols[i].segyData.setIsSampleFormatOverride(b);
        wizard.analyzeHeaders();
    }
*/

    public void setOverrideNSamples(int nSamples)
    {
        overrideNSamples = nSamples;
        segyFormat.setNSamp(nSamples);
        StsSeismicBoundingBox[] vols = wizard.getSelectedVolumes();
    	for(int i=0; i<vols.length; i++)
    		vols[i].segyData.setOverrideNSamples(nSamples);
        wizard.analyzeHeaders();
    }

    public int getOverrideNSamples() { return overrideNSamples; }

    public void setOverrideSampleSpacing(float sampleSpacing)
    {
        overrideSampleSpacing = sampleSpacing;
        segyFormat.setSRate(sampleSpacing);
        StsSeismicBoundingBox[] vols = wizard.getSelectedVolumes();
    	for(int i=0; i<vols.length; i++)
    		vols[i].segyData.setOverrideSampleSpacing(sampleSpacing);
        wizard.analyzeHeaders();
    }

    public float getOverrideSampleSpacing() { return overrideSampleSpacing; }
    
    public void setOverrideHeader(boolean b)
    {
        overrideHeader = b;
        segyFormat.setOverrideHeader(overrideHeader);
        StsSeismicBoundingBox[] vols = wizard.getSelectedVolumes();
    	for(int i=0; i<vols.length; i++)
    		vols[i].segyData.setOverrideHeader(b);
        wizard.analyzeHeaders();
    }

    public boolean getOverrideHeader() { return overrideHeader; }
    
    public void setStartZ(float z)
    {
        startZ = z;
        segyFormat.setStartZ(startZ);
        StsSeismicBoundingBox[] vols = wizard.getSelectedVolumes();
    	for(int i=0; i<vols.length; i++)
    		vols[i].segyData.setStartZ(z);	
    }
    
    public float getStartZ() { return startZ; }
    
    public int getTraceHeaderSize()
    {
        return traceHeaderSize;
    }

    public void setTraceHeaderSize(int traceHeaderSize)
    {
        this.traceHeaderSize = traceHeaderSize;
        segyFormat.setTraceHeaderSize(traceHeaderSize);
        StsSeismicBoundingBox[] vols = wizard.getSelectedVolumes();
    	for(int i=0; i<vols.length; i++)
    		vols[i].setTraceHeaderSize(traceHeaderSize);
    }
    
    public void setTextHeaderFormatString(String string)
    {
        textHeaderFormatString = string;
        segyFormat.setTextHeaderFormatString(string);
        StsSeismicBoundingBox[] vols = wizard.getSelectedVolumes();
    	for(int i=0; i<vols.length; i++)
    		vols[i].segyData.setTextHeaderFormatString(string);
    }
    
    public String getTextHeaderFormatString() { return textHeaderFormatString; }
	
    public int getTextHeaderSize()
    {
        return textHeaderSize;
    }

    public void setTextHeaderSize(int textHeaderSize)
    {
        this.textHeaderSize = textHeaderSize;
        segyFormat.setTextHeaderSize(textHeaderSize);
        StsSeismicBoundingBox[] vols = wizard.getSelectedVolumes();
    	for(int i=0; i<vols.length; i++)
    		vols[i].setTextHeaderSize(textHeaderSize);
    }	
    
    public int getBinaryHeaderSize()
    {
        return binaryHeaderSize;
    }

    public void setBinaryHeaderSize(int binaryHeaderSize)
    {
        this.binaryHeaderSize = binaryHeaderSize;
        segyFormat.setBinaryHeaderSize(binaryHeaderSize);
        StsSeismicBoundingBox[] vols = wizard.getSelectedVolumes();
    	for(int i=0; i<vols.length; i++)
    		vols[i].setBinaryHeaderSize(binaryHeaderSize);
    }
    
    public boolean getIsLittleEndian()
    {
        return isLittleEndian;
    }

    public void setIsLittleEndian(boolean isLittleEndian)
    {
        this.isLittleEndian = isLittleEndian;
        segyFormat.setIsLittleEndian(isLittleEndian);
        StsSeismicBoundingBox[] vols = wizard.getSelectedVolumes();
    	for(int i=0; i<vols.length; i++)
            vols[i].setIsLittleEndian(isLittleEndian);
        wizard.analyzeHeaders();
    }	
    
    public String getOverrideSampleFormatString()
    {
        return overrideSampleFormatString;
    }

    public void setOverrideSampleFormatString(String overrideSampleFormatString)
    {
        this.overrideSampleFormatString = overrideSampleFormatString;
        byte overrideSampleFormat = (byte)StsParameters.getStringMatchIndex(StsSEGYFormat.sampleFormatStrings, overrideSampleFormatString);
        segyFormat.setOverrideSampleFormat(overrideSampleFormat);
        StsSeismicBoundingBox[] vols = wizard.getSelectedVolumes();
    	for(int i=0; i<vols.length; i++)
            vols[i].setOverrideSampleFormat(overrideSampleFormat);
        wizard.analyzeHeaders();
    }

    public String getHeaderSampleFormatString()
    {
        if(wizard == null) return StsSEGYFormat.NONE_STRING;
        StsSeismicBoundingBox[] vols = wizard.getSelectedVolumes();
        if(vols == null || vols.length == 0) return StsSEGYFormat.NONE_STRING;
        int sampleFormat = vols[0].getHeaderSampleFormat();
        for(int i=0; i<vols.length; i++)
        {
            int otherSampleFormat = vols[i].getSampleFormat();
            if(otherSampleFormat != sampleFormat) return StsSEGYFormat.NONE_STRING;
        }
        return StsSEGYFormat.sampleFormatStrings[sampleFormat];
    }

    public boolean getOverride() { return getOverrideHeader(); }
    public void setOverride(boolean b)
    {
        setOverrideHeader(b);
        overrideNSamplesBean.setEditable(b);
        overrideSampleSpacingBean.setEditable(b);

    }

	public void setMessage(String msg)
	{
		messageTxt.setText(msg);
		repaint();
	}

	public void viewTextHeaders()
	{
        if(wizard == null) return;
        Frame frame = wizard.getModel().win3d;

        StsSeismicBoundingBox[] volumes = wizard.getSelectedVolumes();
        if(volumes == null || volumes.length == 0)
        {
            new StsMessage(frame, StsMessage.INFO,  "No volumes are available.");
            return;
        };

        StsTabbedTextAreaDialog dialog = new StsTabbedTextAreaDialog(frame, "SegY Text Headers", false);

		for (int i = 0; i < volumes.length; i++)
		{
			StsTextAreaScrollPane textArea = dialog.addTab(volumes[i].getName());
            String encoder = volumes[i].getSegyFormat().getTextHeaderFormat();
		    String header = volumes[i].segyData.readTextHdr(encoder);
			for (int j = 0; j <= header.length() - LINE_LENGTH; j += LINE_LENGTH)
			{
				String line = header.substring(j, j + LINE_LENGTH-1) + StsParameters.EOL;
				textArea.append(line);
			}
		}

		dialog.setVisible(true);
		dialog.setTabActive(0);
	}

	public void viewBinaryHeaders()
	{
        if(wizard == null) return;
        Frame frame = wizard.getModel().win3d;

        StsSeismicBoundingBox[] volumes = wizard.getSelectedVolumes();
        if(volumes == null || volumes.length == 0)
        {
            new StsMessage(frame, StsMessage.INFO,  "No volumes are available.");
            return;
        };

        StsTabbedTextAreaDialog dialog = new StsTabbedTextAreaDialog(frame, "SegY Binary Headers", false);
		if (volumes == null)
		{
			new StsMessage(frame, StsMessage.WARNING, "No headers to display");
			return;
		}

		for (int i = 0; i < volumes.length; i++)
		{
			StsTextAreaScrollPane textAreaScrollPane = dialog.addTab(volumes[i].getName());
            StsSEGYFormatRec[] allBinaryRecs = volumes[i].getSegyFormat().getAllBinaryRecords();

			for (int n = 0; n < allBinaryRecs.length; n++)
			{
				StsSEGYFormatRec rec = allBinaryRecs[n];
                textAreaScrollPane.append(buildLine(volumes[i].getBinaryHeaderValue(rec), rec.getDescription()) + StsParameters.EOL);
			}
		}

		dialog.setVisible(true);
		dialog.setTabActive(0);
	}

	private String buildLine(double num, String desc)
	{
		StringBuffer lbuf = new StringBuffer(
			"                                                                               ");

		String numText = fmt.format(num);
		if(numText.length() > 15) numText = numText.substring(0, 14);
		lbuf.insert(15 - numText.length(), numText);
		lbuf.insert(17, "-");
		lbuf.insert(19, desc);

		return lbuf.toString();
	}

    static public void main(String[] args)
	{
		com.Sts.MVC.StsModel model = new com.Sts.MVC.StsModel();
		com.Sts.MVC.StsProject project = new com.Sts.MVC.StsProject();
		model.setProject(project);
		StsFileFormatPanel panel = new StsFileFormatPanel(null, true);
		StsJPanel backPanel = StsJPanel.addInsets();
		backPanel.add(panel);
		StsToolkit.createDialog(backPanel, true);
	}
}
