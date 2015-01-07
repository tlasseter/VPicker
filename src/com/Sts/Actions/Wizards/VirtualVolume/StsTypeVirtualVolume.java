package com.Sts.Actions.Wizards.VirtualVolume;

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

public class StsTypeVirtualVolume extends StsWizardStep
{
    StsTypeVirtualVolumePanel panel;
    StsHeaderPanel header;

    public StsTypeVirtualVolume(StsWizard wizard)
	{
        super(wizard);
        panel = new StsTypeVirtualVolumePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Virtual Volume Definition");
        header.setSubtitle("Select Volume Type");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#VirtualVolume");
        header.setInfoText(wizardDialog,"(1) Specify the name of the filter volume\n" +
                           " (2) Select the type of virtual volume that will be defined.\n" +
                           "(3) Press the Next > Button to proceed to the definition screen.\n");
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

    public int getVolumeType()
    {
        return panel.getVolumeType();
    }

    public String getVolumeName()
    {
        return panel.getName();
    }

}
