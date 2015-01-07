package com.Sts.Actions.Wizards.HandVelocity;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

public class StsHandVelocitySurveyDefinition extends StsWizardStep
{
    public StsHandVelocitySurveyDefinitionPanel panel;
    private StsHeaderPanel header;

    public StsHandVelocitySurveyDefinition(StsWizard wizard)
    {
        super(wizard);
        panel = new StsHandVelocitySurveyDefinitionPanel((StsHandVelocityWizard)wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Hand Velocity Survey Definition");
        header.setSubtitle("Define Hand Velocity Geometry");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#HandVel");                        
        header.setInfoText(wizardDialog,"      **** At this point, hand veloicty lines listed in the table need to be associated *****\n" +
                                        "      **** with existing survey lines which contain geometry information for each cdp. *****\n" +
                                        "(1) Select a hand velocity line from the table.\n" +
                                        "(2) Select a corresponding survey definition line containing the cdp geometry information.\n" +
                                        "(2) Verify ALL values in file table are correct before proceeding to processing\n" +
                                        "(3) Press the Finish> Button to complete the association and proceed to Velocity Analysis.\n");

    }

    public boolean start()
    {
        wizard.dialog.setTitle("Define Hand Velocity Geometry");
        panel.initialize();
        return true;
    }

    public boolean end()
    {
/*
        StsSegyLine2d[] lines = null;
		lines = ((StsSegyLine2dWizard)wizard).getSegyVolumesToProcess();
        if( lines == null) return true;

        for(int i = 0; i < lines.length; i++)
		{
            if( lines[i].getOverrideGeom())
            {
                if( !lines[i].checkGeomAttributes())
                    return false;
            }
        }
*/
		return true;
    }

	public void updatePanel()
	{
		panel.updatePanel();
	}
}
