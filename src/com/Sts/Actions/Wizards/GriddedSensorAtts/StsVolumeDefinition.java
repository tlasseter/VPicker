
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.GriddedSensorAtts;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.PreStack3d.StsPreStackSurveyDefinitionPanel;
import com.Sts.Actions.Wizards.PreStack3d.StsPreStackWizard;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.Types.PreStack.*;
import com.Sts.Types.*;

public class StsVolumeDefinition extends StsWizardStep
{
    public StsVolumeDefinitionPanel panel;
    private StsHeaderPanel header;

    public StsVolumeDefinition(StsWizard wizard)
    {
        super(wizard);
        panel = new StsVolumeDefinitionPanel((StsGriddedSensorAttsWizard)wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Volume Definition");
        header.setSubtitle("Define Volume Bounds");
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
        wizard.dialog.setTitle("Define Volume Bounds");
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        ((StsGriddedSensorAttsWizard)wizard).setBoundingBox(panel.getSeismicBoundingBox());
		return true;
    }

}