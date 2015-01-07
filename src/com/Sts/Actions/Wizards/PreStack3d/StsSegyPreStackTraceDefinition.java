package com.Sts.Actions.Wizards.PreStack3d;

import com.Sts.Actions.Wizards.PreStack2d.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.MVC.*;
import com.Sts.Types.PreStack.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

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
abstract public class StsSegyPreStackTraceDefinition extends StsWizardStep
{
	protected StsSegyTraceDefinitionPanel panel;
	protected StsHeaderPanel header;
	protected StsPreStackWizard wizard = null;

    abstract public void analyzeGeometry();

    public StsSegyPreStackTraceDefinition(StsWizard wizard)
	{
		super(wizard);
		this.wizard = (StsPreStackWizard)wizard;

        if( wizard instanceof StsPreStack2dWizard)
            panel = new StsSegyTraceDefinitionPanel2d( wizard, this);
        else
            panel = new StsSegyTraceDefinitionPanel( wizard, this);

        header = new StsHeaderPanel();
		setPanels(panel, header);

		header.setTitle(getTitle());
		header.setSubtitle(getSubtitle());
		header.setLogo(getLogo());
		header.setLink(getLink());
		header.setInfoText(wizardDialog, getInfoText().toString());
	}

	public StsSegyTraceDefinitionPanel getPanel()
	{
		return panel;
	}


	public boolean start()
	{
        panel.initialize();
        panel.updatePanel();
        analyzeGeometry();

        return true;
	}

	public void updatePanel()
	{
		panel.updatePanel();
	}
/**
	 * getAllTraceRecords
	 *
	 * @param volume StsRotatedGridBoundingBox
	 * @return StsSEGYFormatRec[]
	 */
	public StsSEGYFormatRec[] getAllTraceRecords(StsRotatedGridBoundingBox volume)
	{
		return ((StsPreStackSegyLine)volume).getSEGYFormat().getAllTraceRecords();
		//return wizard.getSegyFormat().getAllTraceRecords();
	}

	public StsSEGYFormatRec[] getRequiredTraceRecords(StsRotatedGridBoundingBox volume)
	{
            return ((StsPreStackSegyLine)volume).getSEGYFormat().getRequiredTraceRecords();
	}

/**
	 * getInfoText
	 *
	 * @return StringBuffer
	 */
	public StringBuffer getInfoText()
	{
		StringBuffer sb = new StringBuffer();
		sb.append("      **** Review the trace definition table to ensure trace values are being read correctly. ****\n");
		sb.append("      **** Minimum required fields are the inline and crossline numbers, ensure they are correct. ****\n");
		sb.append("      **** Trace header attributes that will be stored on processing will have grey backgrounds. ****\n");
		sb.append("(1) Set the trace header size if non-standard (240 bytes).\n");
		sb.append("(2) To limit the table to only header attributes that are to be stored press the Selected Attributes box.\n");
		sb.append("(3) Select rows of attributes that require changes and set new format, byte locations and scalar application.\n");
		sb.append("      ***** Scalar is in location 69 and is multiplied if positive, divided if negative per SEGY standard ****\n");
		sb.append("(4) Press Update button to save the edits to the selected trace header attribute definition.\n");
		sb.append("      ***** Any changes will result in a re-scan of the headers to see if attribute is being correctly read. ****\n");
		sb.append("(5) Add additional desired attributes by selecting the row and pressing the Add button.\n");
		sb.append("(6) Remove attributes by changing table to selected attributes, selecting attribute and pressing remove button.\n");
		sb.append("      ***** CHANGES WILL APPLY TO ALL THE FILES IN THE TABLE AT THE BOTTOM ON THE SCREEN. ****\n");
		sb.append("(7) Once the trace header information has been reviewed and values set, press the Next>> Button.\n");
		sb.append("      ***** All values in file table should be correct except coordinates prior to moving to next screen. ****\n");
		return sb;
	}

/**
	 * getLink
	 *
	 * @return String
	 */
	public String getLink()
	{
		return "http://www.s2ssystems.com/marketing/s2ssystems/SegyLoad.html";
	}

/**
	 * getLogo
	 *
	 * @return String
	 */
	public String getLogo()
	{
		return "TsunamiLogo.gif";
	}

/**
	 * getNTraces
	 *
	 * @param volume StsRotatedGridBoundingBox
	 * @return int
	 */
	public int getNTraces(StsRotatedGridBoundingBox volume)
	{
		return ((StsPreStackSegyLine)volume).getTotalTraces();
	}

/**
	 * getOverrideHeaderAllowed
	 *
	 * @return boolean
	 */
	public boolean getOverrideHeaderAllowed()
	{
		return false;
	}

/**
	 * getSubtitle
	 *
	 * @return String
	 */
	public String getSubtitle()
	{
		return "Trace Definition";
	}

/**
	 * getTitle
	 *
	 * @return String
	 */
	public String getTitle()
	{
		return "Pre-Stack3d Details";
	}

/**
	 * getTraceHeaderBinary
	 *
	 * @param volume StsRotatedGridBoundingBox
	 * @param nTrace int
	 * @return byte[]
	 */
	public byte[] getTraceHeaderBinary(StsRotatedGridBoundingBox volume, int nTrace)
	{
		return ((StsPreStackSegyLine)volume).getTraceHeaderBinary(nTrace);
	}

/**
	 * getTraceHeaderSize
	 *
	 * @param volume StsRotatedGridBoundingBox
	 * @return long
	 */
	public long getTraceHeaderSize(StsRotatedGridBoundingBox volume)
	{
		return ((StsPreStackSegyLine)volume).getSEGYFormat().getTraceHeaderSize();
		//return wizard.getSegyFormat().getTraceHeaderSize();
	}

/**
	 * overrides super.end()
	 */
	public boolean end()
	{
		return true;
	}

	// overrides super.initSetValues(StsRotatedGridBoundingBox volume)
	public void initSetValues(StsRotatedGridBoundingBox volume)
	{
		((StsPreStackSegyLine)volume).initializeSegyIO();
	}

	public boolean hasGeometry()
	{
		return true;
	}

	public int[] getTraceStartAndCount(StsRotatedGridBoundingBox volume, int traceNumber)
	{
		int nTracesPerGather = ((StsPreStackSegyLine)volume).getNTracesPerGather();
		if (((StsPreStackSegyLine)volume).isTracesPerGatherSame)
		{
			traceNumber = StsMath.intervalRoundDown(traceNumber, nTracesPerGather);
		}
		else
		{
			int[] gatherTraceRange = ((StsPreStackSegyLine)volume).findGatherTraceRange(traceNumber, nTracesPerGather);
			if (gatherTraceRange == null)
				return new int[]
					{0, 50};
			traceNumber = gatherTraceRange[0];
			int nLastTrace = gatherTraceRange[1];
			nTracesPerGather = nLastTrace - traceNumber + 1;
		}
		return new int[]
			{traceNumber, nTracesPerGather};
	}

	public boolean getOverrideStep()
	{
		return false;
	}

	// JKF -- not sure if this method is required...evaluate with Tom
	public boolean saveTemplate(StsSegyVolume vol)
	{
		String templateName = null;

		try
		{
			StsTextAreaDialog dialog = new StsTextAreaDialog(model.win3d, "Enter Template Name...", null, 1, 40);
			dialog.setVisible(true);
			templateName = dialog.getText();
			if ((templateName == null) || (templateName.length() < 1))
			{
				StsMessageFiles.infoMessage("Invalid Template Name");
				return false;
			}
			templateName = StsStringUtils.cleanString(templateName);
			String pathname = model.getProject().getDataFullDirString() + "seis3dTemplate.txt." + templateName;
			StsParameterFile.writeObjectFields(pathname, vol, StsSegyVolume.class, StsBoundingBox.class);
			return true;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSegyTraceDefinition.saveTemplate() failed.",
										 e, StsException.WARNING);
			return false;
		}
	}

	private StsRotatedGridBoundingBox[] getSegyVolumes()
	{
		return (wizard).getVolumes();
	}

}
