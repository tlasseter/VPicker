package com.Sts.Actions.Wizards.FractureTrack;

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

public class StsDefineTracker extends StsWizardStep
{
	StsDefineTrackerPanel panel;
    StsHeaderPanel header;

    public StsDefineTracker(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineTrackerPanel((StsFractureTrackWizard)wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Fracture Tracking");
        header.setSubtitle("Define Tracking Criteria");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#FractureTrack");        
        header.setInfoText(wizardDialog,"(1) \n" +
        		"(5) Press the Finish Button, when done.\n");
        header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/FractureTrack.html");
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
