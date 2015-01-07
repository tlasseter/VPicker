package com.Sts.Actions.Wizards.PreStack3d;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.Types.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2007</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

public class StsPreStackBatch extends StsWizardStep
{
	private StsPreStackBatchPanel panel;
	private StsHeaderPanel header;

	public StsPreStackBatch(StsWizard wizard)
	{
		super(wizard);
		panel = new StsPreStackBatchPanel(wizard, this);
		header = new StsHeaderPanel();
		setPanels(panel, header);
		header.setTitle("Pre-Stack Volumes Set to Process");
		header.setSubtitle("Add Volumes to Batch Process List");
		header.setInfoText(wizardDialog,"(1) Add files to the batch process list, if desired.\n" +
										"      **** All Files in the top list will be add in the bottom list *****\n" +
										"(2) Remove files from the bottom, awaiting processing, list if they are not to be processed.\n" +
										"(3) Use the Back> Button to go back to the start and add more files for processing.\n" +
										"(4) Press the Next> Button to proceed to processing the files.\n");
	}

	public StsPreStackBatchPanel getPanel()
	{
        wizard.dialog.setTitle("Crop Optimized Volumes");
		return panel;
	}

	public boolean start()
	{
		panel.initialize();
		return true;
	}

	public StsSeismicBoundingBox[] getSegyVolumes()
	{
		return ((StsPreStackWizard)wizard).getSegyVolumes();
	}

	public boolean end()
	{
		return true;
	}

	public void updatePanel()
	{
		panel.updatePanel();
	}

	public void addVolumesForProcessing()
	{
		((StsPreStackWizard)wizard).addSegyVolumesToProcessList();
		updatePanel();
	}
}
