package com.Sts.Actions.Wizards.CombinationVolume;

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

public class StsDefineCombinationVolume extends StsWizardStep
{
    StsDefineCombinationVolumePanel panel;
    StsHeaderPanel header;

    public StsDefineCombinationVolume(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineCombinationVolumePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 500));
        header.setTitle("PostStack3d Combination Definition");
        header.setSubtitle("Defining PostStack3d Type");
        header.setInfoText(wizardDialog,"(1) Specify the name of the combination volume being produced.\n" +
                           "(2) Select the type of volume to produce.\n" +
                           "     **** A description of each type is provided upon selection of the type. ****\n" +
                           "     **** Pointsets can only be created using the multi-volume math option. ****\n" +
                           "(3) If the selected type is box set summation, specify whether the resulting\n" +
                           "    volume includes data inside boxes or outside using the Inclusive checkbox.\n" +
                           "(4) Press the Next>> Button.\n");
        header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/VolumeCombine.html");
    }

    public boolean start()
    {
        return true;
    }

    public boolean end()
    {
        return true;
    }

    public boolean isInclusive() { return panel.isInclusive(); }
    public byte getType() { return panel.getType(); }
    public String getVolumeName() { return panel.getVolumeName(); }
}
