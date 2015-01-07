package com.Sts.Actions.Wizards.StsObjectLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;

public class StsObjectSelect extends StsWizardStep {
	StsObjectSelectPanel panel;
	StsHeaderPanel header;

	public StsObjectSelect(StsObjectLoadWizard wizard)
	{
		super(wizard);

		panel = new StsObjectSelectPanel(wizard, this);
		header = new StsHeaderPanel();
		setPanels(panel, header);
		header.setTitle("StsObject Load");
		header.setSubtitle("Select StsObject Files");
        panel.setPreferredSize(new Dimension(500, 300));
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#StsObjectLoad");
        header.setInfoText(wizardDialog,"(1) Navigate to the directory containing the StsObject Files using the Dir Button.\n" +
                                        "(2) Select files from the left list and place them in the right using the controls between the lists.\n" +
                                        "(3) Press the Next> Button");
	}

	public StsObjectSelectPanel getPanel()
	{
		return panel;
	}

	public boolean start()
	{
		StsObjectLoadWizard sensorWizard = (StsObjectLoadWizard)wizard;
		panel.initialize();
		return true;
	}

	public boolean end()
	{
		try
		{
			if (((StsObjectLoadWizard)wizard).getObjectFiles() == null)
			{
				new StsMessage(wizard.frame, StsMessage.ERROR,
							   "No sensor files selected: select or cancel.");
				return false;
			}
			return true;
		}
		catch (Exception e)
		{
			StsException.outputException("StsObjectSelect.end() failed.", e, StsException.WARNING);
			return false;
		}
	}
}
