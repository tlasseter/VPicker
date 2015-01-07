package com.Sts.Actions.Wizards.VirtualVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsDefineFloatVolume extends StsWizardStep
{
	StsDefineFloatVolumePanel panel;
    StsHeaderPanel header;
	
    public StsDefineFloatVolume(StsWizard wizard)
	{
        super(wizard);
        panel = new StsDefineFloatVolumePanel((StsVirtualVolumeWizard)wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 200));
        header.setTitle("Virtual Volume Definition");
        header.setSubtitle("Define Float Resolution Parameters");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#VirtualVolume");
        header.setInfoText(wizardDialog,"(1) Specify the value that you want to be considered a null.\n" +
                           "   **** Example: Difference two volumes and specify zero as null value\n" +
                           "   **** will result in all zeroes being transparent.\n" +
                           "(3) Press the Finish> Button to create the volume.\n");
    }

    public boolean start()
	{
        return true;
    }

    public boolean end()
    {
        return true;
    }


}
