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

public class StsTypeFilterVolume extends StsWizardStep
{
    StsTypeFilterVolumePanel panel;
    StsHeaderPanel header;

    public StsTypeFilterVolume(StsWizard wizard)
	{
        super(wizard);
        panel = new StsTypeFilterVolumePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 200));
        header.setTitle("Filter PostStack3d Definition");
        header.setSubtitle("Select Filter PostStack3d Type");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#VirtualVolume");
        header.setInfoText(wizardDialog,"(1) Specify the name of the virtual volume\n" +
                           " (2) Select the type of virtual volume that will be defined.\n" +
                           "(3) Press the Next > Button to proceed to the definition screen.\n");
    }

    public boolean start()
    {
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
