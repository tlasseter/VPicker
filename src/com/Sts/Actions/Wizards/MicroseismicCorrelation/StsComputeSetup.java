package com.Sts.Actions.Wizards.MicroseismicCorrelation;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;

import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsComputeSetup extends StsWizardStep
{
    StsComputeSetupPanel panel;
    StsHeaderPanel header;
	StsSeismicVolume[] volumes = null;

    public StsComputeSetup(StsWizard wizard)
	{
        super(wizard);
        panel = new StsComputeSetupPanel((StsMicroseismicCorrelationWizard)wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 200));
        header.setTitle("Compute Correlation Attributes");
        header.setSubtitle("Defining Correlation Parameters");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#MicroseismicCorrelation");                                
        header.setInfoText(wizardDialog,"(1) Select the volume to operate on.\n" +
                           "(2) Define the operator dimensions.\n" +
                           "(3) Press the Finish> Button to compute attribute\n");
    }

    public boolean start()
	{
        panel.initialize();
        return true;
    }

    public boolean end()
    {
    	return true;
    }
}