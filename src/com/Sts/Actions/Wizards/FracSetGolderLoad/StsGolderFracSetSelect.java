package com.Sts.Actions.Wizards.FracSetGolderLoad;

import com.Sts.Actions.Wizards.FracSetGolderLoad.StsGolderFracSetLoadWizard;
import com.Sts.Actions.Wizards.FracSetGolderLoad.StsGolderFracSetSelectPanel;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.Actions.Wizards.WizardHeaders.StsHeaderPanel;
import com.Sts.UI.StsMessage;
import com.Sts.Utilities.StsException;

import java.awt.*;

public class StsGolderFracSetSelect extends StsWizardStep
{
	StsGolderFracSetSelectPanel panel;
	StsHeaderPanel header;

	public StsGolderFracSetSelect(StsGolderFracSetLoadWizard wizard)
	{
		super(wizard);

		panel = new StsGolderFracSetSelectPanel(wizard, this);
		header = new StsHeaderPanel();
		setPanels(panel, header);
		header.setTitle("FracSet Load");
		header.setSubtitle("Select FracSet Files");
        panel.setPreferredSize(new Dimension(500, 300));
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#FracSet");
        header.setInfoText(wizardDialog,"(1) Navigate to the directory containing the FracSet Files using the Dir Button.\n" +
                                        "      **** All csv and txt files in the selected directory will be placed in the left list.\n" +
                                        "(2) Select files from the left list and place them in the right using the controls between the lists.\n" +
                                        "(3) Press the Next> Button");
	}

	public StsGolderFracSetSelectPanel getPanel()
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
			if (((StsGolderFracSetLoadWizard)wizard).getFracSetFiles() == null)
			{
				new StsMessage(wizard.frame, StsMessage.ERROR,
							   "No FracSet files selected: select or cancel.");
				return false;
			}
			return true;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSFracSetSelect.end() failed.", e, StsException.WARNING);
			return false;
		}
	}
}
