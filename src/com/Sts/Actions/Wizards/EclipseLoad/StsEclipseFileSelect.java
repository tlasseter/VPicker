package com.Sts.Actions.Wizards.EclipseLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;

public class StsEclipseFileSelect extends StsWizardStep
{
	StsEclipseFileSelectPanelNew panel;
	StsHeaderPanel header;

	public StsEclipseFileSelect(StsEclipseLoadWizard wizard)
	{
		super(wizard);

		panel = new StsEclipseFileSelectPanelNew(wizard, this);
		header = new StsHeaderPanel();
		setPanels(panel, header);
		header.setTitle("Eclipse Restart File Load");
		header.setSubtitle("Select Restart File");
        panel.setPreferredSize(new Dimension(500, 300));
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#EclipseLoad");
        header.setInfoText(wizardDialog,"(1) Navigate to the directory containing the Eclipse Restart Files using the Dir Button.\n" +
                                        "      **** All unrst files in the selected directory will be placed in the left list.\n" +
                                        "(2) Select files from the left list and place them in the right using the controls between the lists.\n" +
                                        "(3) Press the Next> Button");
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
			if (((StsEclipseLoadWizard)wizard).getRestartFile() == null)
			{
				new StsMessage(wizard.frame, StsMessage.ERROR,  "No Eclipse restart files selected: select or cancel.");
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