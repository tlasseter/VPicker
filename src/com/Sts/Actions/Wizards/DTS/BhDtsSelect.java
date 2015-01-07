package com.Sts.Actions.Wizards.DTS;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;

public class BhDtsSelect extends StsWizardStep
{
	BhDtsSelectPanel panel;
	StsHeaderPanel header;

	public BhDtsSelect(BhDtsLoadWizard wizard)
	{
		super(wizard);

		panel = new BhDtsSelectPanel(wizard, this);
		header = new StsHeaderPanel();
		setPanels(panel, header);
		header.setTitle("Distributed Temperature Sensor Load");
		header.setSubtitle("Select DTS Files");
        panel.setPreferredSize(new Dimension(500, 300));
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#DTS");
        header.setInfoText(wizardDialog,"(1) Navigate to the directory containing the Sensor Files using the Dir Button.\n" +
                                        "      **** All csv and txt files in the selected directory will be placed in the left list.\n" +
                                        "(2) Select files from the left list and place them in the right using the controls between the lists.\n" +
                                        "(3) Press the Next> Button");
	}

	public BhDtsSelectPanel getPanel()
	{
		return panel;
	}

	public boolean start()
	{
		panel.initialize();
		return true;
	}

	public boolean end()
	{
		try
		{
			if (((BhDtsLoadWizard)wizard).getFiles() == null)
			{
				new StsMessage(wizard.frame, StsMessage.ERROR,  "No DTS files selected: select or cancel.");
				return false;
			}
			return true;
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "end", e);
			return false;
		}
	}
}
