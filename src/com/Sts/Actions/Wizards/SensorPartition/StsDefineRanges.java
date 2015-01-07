package com.Sts.Actions.Wizards.SensorPartition;

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.Actions.Wizards.WizardHeaders.StsHeaderPanel;

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
        panel = new StsDefineRangesPanel((StsSensorPartitionWizard)wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Analyze Ranges");
        header.setSubtitle("Define Criteria & Analyze Results");
        header.setInfoText(wizardDialog,"(1) Set the min and max time values for each partition.\n" +
        		"(3) Press the Finish Button to create the newly partitioned sensors and exit the wizard.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#SensorPartition");
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