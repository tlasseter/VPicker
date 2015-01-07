package com.Sts.Actions.Wizards.ProximityAnalysis;

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

public class StsDefineProximity extends StsWizardStep
{
	StsDefineProximityPanel panel;
    StsHeaderPanel header;

    public StsDefineProximity(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineProximityPanel((StsProximityAnalysisWizard)wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Analyze Proximity");
        header.setSubtitle("Define Criteria & Analyze Results");
        header.setInfoText(wizardDialog,"(1) Specify the maximum distance from a well.\n" +
        		"(2) Specify whether distances are computed in 2D or 3D space.\n" +
        		"(3) Press the Run Analysis Button to compute and display results.\n" +
        		"(4) Press the Export View Button to export the events in view to a file.\n" +
        		"(5) Press the Finish Button, when done.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#ProximityAnalysis");
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
