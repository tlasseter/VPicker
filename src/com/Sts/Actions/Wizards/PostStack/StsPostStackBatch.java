package com.Sts.Actions.Wizards.PostStack;

import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

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

public class StsPostStackBatch extends StsWizardStep
{
	protected StsPostStackBatchPanel panel;
	protected StsHeaderPanel header;

	public StsPostStackBatch(StsWizard wizard)
	{
		super(wizard);
		panel = new StsPostStackBatchPanel(wizard, this);
		header = new StsHeaderPanel();
		setPanels(panel, header);
		header.setTitle("Post-stack Volumes Set to Process");
		header.setSubtitle("Add Volumes to Batch Process List");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#ProcessSeismic");
		header.setInfoText(wizardDialog,"(1) Add files to the batch process list, if desired.\n" +
										"      **** All Files in the top list will be add in the bottom list *****\n" +
										"(2) Remove files from the bottom, awaiting processing, list if they are not to be processed.\n" +
										"(3) Use the Back> Button to go back to the start and add more files for processing.\n" +
										"(4) Press the Next> Button to proceed to processing the files.\n");
	}

	public StsPostStackBatchPanel getPanel()
	{
        wizard.dialog.setTitle("Crop Optimized Volumes");
		return panel;
	}

	public boolean start()
	{
		panel.initialize();
		return true;
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
		((StsSeismicWizard)wizard).addSegyVolumesToProcessList();
		updatePanel();
	}
}
