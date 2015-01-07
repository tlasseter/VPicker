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

public class StsDefineDifferenceVolume extends StsWizardStep
{
    StsDefineDifferenceVolumePanel panel;
    StsHeaderPanel header;

    public StsDefineDifferenceVolume(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineDifferenceVolumePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 500));
        header.setTitle("PostStack3d Combination Definition");
        header.setSubtitle("Defining Individual Difference PostStack3d");
        header.setInfoText(wizardDialog,"(1) Select the two volumes to be differenced.\n" +
                                   "(2) Specify the data range to be used for scaling the combination volume.\n" +
                                   "     **** PostStack3d one data histogram is isVisible, switch volumes to view the\n" +
                                   "          data distribution of the second volume.\n" +
                                   "(3) Press the Next>> Button.\n");
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

    public StsSeismicVolume getVolumeOne() { return panel.getSelectedSeismicOneVolume(); }
    public StsSeismicVolume getVolumeTwo() { return panel.getSelectedSeismicTwoVolume(); }
    public float getDataMin() { return panel.getDataMin(); }
    public float getDataMax() { return panel.getDataMax(); }

}
