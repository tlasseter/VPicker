package com.Sts.Actions.Wizards.PostStack3d;

import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

public class StsPostStackSurveyDefinition3d extends StsWizardStep
{
    public StsPostStackSurveyDefinitionPanel3d panel;
    private StsHeaderPanel header;

    public StsPostStackSurveyDefinition3d(StsWizard wizard)
    {
        super(wizard);
        panel = new StsPostStackSurveyDefinitionPanel3d((StsSeismicWizard)wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Post-Stack3d Survey Definition");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#PostStack3d");                                
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
		return true;
    }
}
