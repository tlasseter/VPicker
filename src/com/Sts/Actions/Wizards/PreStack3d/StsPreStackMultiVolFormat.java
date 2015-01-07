package com.Sts.Actions.Wizards.PreStack3d;

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

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

public class StsPreStackMultiVolFormat extends StsWizardStep
{
	private StsPreStackMultiVolFormatPanel panel;
	private StsHeaderPanel header;

	public StsPreStackMultiVolFormat(StsWizard wizard)
	{
		super(wizard);
		panel = new StsPreStackMultiVolFormatPanel(wizard, this);
		header = new StsHeaderPanel();
		setPanels(panel, header);
		header.setTitle("Pre-Stack Multi-PreStack3d Definition");
		header.setSubtitle("Verify the Multi-PreStack3d Definition");
		header.setInfoText(wizardDialog,"(1) ");
	}

	public StsPreStackMultiVolFormatPanel getPanel()
	{
		return panel;
	}

	public boolean start()
	{
        wizard.dialog.setTitle("Verify the Multi-PreStack3d Definition");
		panel.initialize();
		return true;
	}

	public boolean end()
	{
		return true;
	}

	public void updatePanel()
	{
	}
}

