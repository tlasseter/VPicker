package com.Sts.Actions.Wizards.CombinationVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;

import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsDefineCummulativeDifferenceVolume extends StsWizardStep
{
    StsDefineCummulativeDifferenceVolumePanel panel;
    StsHeaderPanel header;

    public StsDefineCummulativeDifferenceVolume(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineCummulativeDifferenceVolumePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 500));
        header.setTitle("PostStack3d Definition");
        header.setSubtitle("Defining Cumulative Difference PostStack3d");
        header.setInfoText(wizardDialog,"(1) Select the desired volumes from the left list and move them to the right list.\n" +
                                   "(2) Position the volumes as desired using the up and down arrows.\n" +
                                   "     **** The first volume will be differenced from the second and summed\n" +
                                   "          with the difference of the third from fourth volume, etc...\n" +
                                   "(3) When a volume is selected, a data histogram will be isVisible allowing\n" +
                                   "    the setting of data scaling for the resulting combination volume.\n" +
                                   "(4) Press the Next>> Button.\n");
        header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/VolumeCombine.html");
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

//    public boolean getAutoChk() { return panel.getAutoChk(); }
    public float getDataMin() { return panel.getDataMin(); }
    public float getDataMax() { return panel.getDataMax(); }
    public StsSeismicVolume[] getVolumes()
    {
        return panel.getVolumeList();
    }
}
