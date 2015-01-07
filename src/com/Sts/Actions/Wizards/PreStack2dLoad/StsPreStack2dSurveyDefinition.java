
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.PreStack2dLoad;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.Types.*;

public class StsPreStack2dSurveyDefinition extends StsWizardStep
{
	public StsPreStack2dSurveyDefinitionPanel panel;
	private StsHeaderPanel header;

	public StsPreStack2dSurveyDefinition(StsWizard wizard)
	{
		super(wizard);
		panel = new StsPreStack2dSurveyDefinitionPanel((StsPreStackLoad2dWizard)wizard, this);
		header = new StsHeaderPanel();
		setPanels(panel, header);
		header.setTitle("Pre-Stack 2D lines Survey Definition");
		header.setSubtitle("Define Survey (Velocity Model) Bounds");
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
        wizard.dialog.setTitle("Define Survey (Velocity Model) Bounds");
		panel.initialize();
		//wizard.enableFinish();
		return true;
	}

	public boolean end()
	{
		//((StsPreStackLoad2dWizard)wizard).getVolumes()[0].addVolumeBoundingBox(panel.getSurveyDefinitionBoundingBox());
        //((StsPreStackLoad2dWizard)wizard).getVolumes()[0].setIsMigrated(panel.getMigrated());
		return true;
	}

	public StsRotatedGridBoundingBox getBoundingBox()
	{
		return panel.getSurveyDefinitionBoundingBox();
    }
	public boolean getOverrideGeometry()
	{
		return panel.getOverride();
	}
}
