package com.Sts.Actions.Wizards.PostStack2d;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.Types.*;

public class StsPostStackSurveyDefinition2d extends StsWizardStep
{
    public StsPostStackSurveyDefinitionPanel2d panel;
    private StsHeaderPanel header;

    public StsPostStackSurveyDefinition2d(StsWizard wizard)
    {
        super(wizard);
        panel = new StsPostStackSurveyDefinitionPanel2d((StsPostStack2dWizard)wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Post-stack Line Survey Definition");
        header.setSubtitle("Define Survey Line Associations");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#PostStack2d");                 
        header.setInfoText(wizardDialog,"      **** At this point, survey lines listed in the table need to be associated with *****\n" +
                                        "      **** existing survey lines which contain geometry information for each cdp. *****\n" +
                                        "(1) Select an input line from the table.\n" +
                                        "(2) Select a corresponding survey definition line containing the cdp geometry information.\n" +
                                        "(3) Verify ALL values in file table are correct before proceeding to processing\n" +
                                        "(4) Press the Next> Button to proceed to verification and processing screen.\n");

    }

    public boolean start()
    {
        wizard.dialog.setTitle("Define Survey Line Associations");
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        StsSeismicBoundingBox[] lines = ((StsPostStack2dWizard)wizard).getSegyVolumesToProcess();
        if( lines == null) return true;

        for(int i = 0; i < lines.length; i++)
		{
            StsSegyLine2d line = (StsSegyLine2d)lines[i];
            if(line.assocLine != null)
            {
                if(!line.checkGeomAttributes())
                    return false;
            }
        }

/*
		if(panel.getUseDecimation())
		{
			for(int i = 0; i < volumes.length; i++)
			{
				int inlineDecimation = panel.getInlineDecimation();
				int xlineDecimation = panel.getXlineDecimation();
				volumes[i].setInlineDecimation(inlineDecimation);
				volumes[i].setXlineDecimation(xlineDecimation);
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
