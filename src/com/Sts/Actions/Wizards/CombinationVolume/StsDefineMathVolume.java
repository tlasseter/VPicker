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

public class StsDefineMathVolume extends StsWizardStep
{
    StsDefineMathVolumePanel panel;
    StsHeaderPanel header;

    public StsDefineMathVolume(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineMathVolumePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 500));
        header.setTitle("PostStack3d Definition");
        header.setSubtitle("Defining Math PostStack3d");
        header.setInfoText(wizardDialog,"(1) Select the desired volumes from the left list and move them to the right list.\n" +
                                   "(2) Position the volumes as desired using the up and down arrows or sort combolist.\n" +
                                   "(3) When a volume is selected, a data histogram will be isVisible allowing\n" +
                                   "    the setting of data scaling for the resulting combination volume.\n" +
                                   "   ***** The keyboard arrow keys can be used to move through the list and see each histogram.\n" +
                                   "(4) If a point set is desired select the checkbox and input a name.\n" +
                                   "   ***** A point set is a group of samples extracted from each volume with \n" +
                                   "         the volume and amplitude. Once built it can be animated to view changes \n" +
                                   "         through time. ***** \n" +
                                   "(5) To reset the histogram to full range, press the Reset Histogram button.\n" +
                                   "(6) Select the operation to be performed on the selected volumes.\n" +
                                   "(7) Press the Next>> Button.\n");
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
    public float getDataMin() { return panel.getDataMin(); }
    public float getDataMax() { return panel.getDataMax(); }

    public StsSeismicVolume[] getVolumes()
    {
        return panel.getVolumeList();
    }
    public byte getOperator()
    {
        return panel.getOperator();
    }

}
