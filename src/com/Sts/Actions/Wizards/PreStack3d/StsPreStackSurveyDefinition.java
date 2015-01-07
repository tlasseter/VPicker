
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.PreStack3d;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.Types.PreStack.*;
import com.Sts.Types.*;

public class StsPreStackSurveyDefinition extends StsWizardStep
{
    public StsPreStackSurveyDefinitionPanel panel;
    private StsHeaderPanel header;

    public StsPreStackSurveyDefinition(StsWizard wizard)
    {
        super(wizard);
        panel = new StsPreStackSurveyDefinitionPanel((StsPreStackWizard)wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Pre-Stack3d Survey Definition");
        header.setSubtitle("Define Survey Bounds");
        header.setInfoText(wizardDialog,"      **** At this point, all information in table at bottom of screen should be correct *****\n" +
                                        "      **** with the possible exception of coordinates which can be overriden on this screen. *****\n" +
                                        "(1) If required, select the override check, and viable override methods and fields will be activated.\n" +
                                        "(2) The user can either manually input values or use the current project definition.\n" +
                                        "      **** The software will decide if project override is possible. If not, option is not shown.\n" +
                                        "(3) If user specified, supply origin, interval and angle based on supplied starting inline, crossline.\n" +
                                        "(4) Verify ALL values in file table are correct before proceeding to processing\n" +
                                        "(5) Press the Next> Button to proceed to verification and processing screen.\n");

    }
    public boolean start()
    {
        wizard.dialog.setTitle("Define Survey Bounds");
        panel.initialize();
        return true;
    }

    public boolean end()
    {
		StsSeismicBoundingBox[] volumes = null;
		volumes = ((StsPreStackWizard)wizard).getSegyVolumesToProcess();
		for(int i = 0; i < volumes.length; i++)
		{
			((StsPreStackSegyLine)volumes[i]).setOverrideGeometry(panel.getSurveyDefinitionBoundingBox());
		}
		if(panel.getUseDecimation())
		{
			for(int i = 0; i < volumes.length; i++)
			{
                StsPreStackSegyLine volume = (StsPreStackSegyLine)volumes[i];
                int inlineDecimation = panel.getInlineDecimation();
				int xlineDecimation = panel.getXlineDecimation();
				volume.setInlineDecimation(inlineDecimation);
				volume.setXlineDecimation(xlineDecimation);
			}
		}
		return true;
    }

	public boolean getOverrideGeometry()
	{
		return panel.getOverride();
	}

	public void updatePanel()
	{
		panel.updatePanel();
	}
}
