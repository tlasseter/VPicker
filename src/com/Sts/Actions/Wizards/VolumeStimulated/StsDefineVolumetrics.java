package com.Sts.Actions.Wizards.VolumeStimulated;

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

public class StsDefineVolumetrics extends StsWizardStep
{
	StsDefineVolumetricsPanel panel;
    StsHeaderPanel header;

    public StsDefineVolumetrics(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineVolumetricsPanel((StsVolumeStimulatedWizard)wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Analyze Volumetrics");
        header.setSubtitle("Define Analysis Criteria");
        header.setInfoText(wizardDialog,"(1) Specify the fracture azimuth.\n" +
        		"(2) Specify the ratio along azimuth versus perpendicular to azimuth.\n" +
        		"(3) Press the Run Analysis Button to compute volume and display results.\n" +
        		"(4) Press the Finish Button, when done.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#VolumeStimulated");
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

}
