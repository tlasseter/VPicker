package com.Sts.Actions.Wizards.SensorLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;

public class StsSensorSelect extends StsWizardStep {
	StsSensorSelectPanel panel;
	StsHeaderPanel header;

	public StsSensorSelect(StsSensorLoadWizard wizard)
	{
		super(wizard);

		panel = new StsSensorSelectPanel(wizard, this);
		header = new StsHeaderPanel();
		setPanels(panel, header);
		header.setTitle("Sensor Load");
		header.setSubtitle("Select Sensor Files");
        panel.setPreferredSize(new Dimension(500, 300));
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#SensorLoad");
        header.setInfoText(wizardDialog,"(1) Navigate to the directory containing the Sensor Files using the Dir Button.\n" +
                                        "      **** All csv and txt files in the selected directory will be placed in the left list.\n" +
                                        "(2) Select files from the left list and place them in the right using the controls between the lists.\n" +
                                        "(3) Press the Next> Button");
	}

	public StsSensorSelectPanel getPanel()
	{
		return panel;
	}

	public boolean start()
	{
		StsSensorLoadWizard sensorWizard = (StsSensorLoadWizard)wizard;
		panel.initialize();
		return true;
	}

	public boolean end()
	{
		try
		{
			if (((StsSensorLoadWizard)wizard).getSensorFiles() == null)
			{
				new StsMessage(wizard.frame, StsMessage.ERROR,
							   "No sensor files selected: select or cancel.");
				return false;
			}
			return true;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSSensorSelect.end() failed.", e, StsException.WARNING);
			return false;
		}
	}
}
