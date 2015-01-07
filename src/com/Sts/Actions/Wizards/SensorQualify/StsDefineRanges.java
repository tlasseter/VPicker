package com.Sts.Actions.Wizards.SensorQualify;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsDefineRanges extends StsWizardStep
{
	StsDefineRangesPanel panel;
    StsHeaderPanel header;

    public StsDefineRanges(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineRangesPanel((StsSensorQualifyWizard)wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Analyze Ranges");
        header.setSubtitle("Define Criteria & Analyze Results");
        header.setInfoText(wizardDialog,"(1) Set the min and max values for the attributes.\n" +
        		" **** Supplied attributes are those common to all selected sensors. ****\n" +
        		"(2) Press the Export View Button to export the events in view to a file.\n" +
        		"(3) Press the Finish Button, when done.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#SensorQualify");
    }

    public boolean start()
    {
        panel.initialize();
        wizard.enableFinish();
        return true;
    }

    public boolean end()
    {
        return true;
    }
    
    public void setAnalysisMessage(String msg)
    {
    	panel.setMessage(msg);
    }
}
