package com.Sts.Actions.Wizards.FracSim;

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

public class StsDefineFracSimVolume extends StsWizardStep
{
    StsDefineFracSimVolumePanel panel;
    StsHeaderPanel header;

    public StsDefineFracSimVolume(StsFracSimWizard wizard)
    {
        super(wizard);
        panel = new StsDefineFracSimVolumePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 500));
        header.setTitle("Sub-PostStack3d Definition");
        header.setSubtitle("Defining Sub-PostStack3d Type");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#FracSim");
        header.setInfoText(wizardDialog,"(1) Select the sub-volume type that is being defined.\n" +
                           "    **** Single Surface is defined with a surface and an offset above and below the surface.  ****\n" +
                           "    **** Dual Surface is defined as the data between two surfaces with optional offsets.  ****\n" +
                           "    **** Box Set is defined graphically by selecting points and setting the size interactively. ****\n" +
                           "    **** Well Set is defined with a well and a starting and ending offset and radius to define a cylinder.  ****\n" +
                           "    **** Reservoir Unit is the data contained in a pre-defined reservoir unit  ****\n" +
                           "(2) Press the Next> Button to proceed to the specific type definition screen.\n");
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