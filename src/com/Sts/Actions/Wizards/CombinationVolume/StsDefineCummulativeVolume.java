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

public class StsDefineCummulativeVolume extends StsWizardStep
{
    StsDefineCummulativeVolumePanel panel;
    StsHeaderPanel header;

    static final int SUM = 0;
    static final int AVERAGE = 1;
    static final int FIRST = 2;
    static final int LAST = 3;

    public StsDefineCummulativeVolume(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineCummulativeVolumePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 500));
        header.setTitle("PostStack3d Definition");
        header.setSubtitle("Defining Cumulative Difference PostStack3d");
        header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/VolumeCombine.html");
        header.setInfoText(wizardDialog,"(1) Select all the subVolumes that you want included in cumulative volume\n" +
                           "(2) Select overlap mode\n" +
                           "    **** Sum mode will add overlapping areas **** \n" +
                           "    **** Average mode will average all overlapping data areas **** \n" +
                           "    **** First will accept the first data encountered and ignore any other ****\n" +
                           "    **** Last will accept the last data encountered and ignore any other ****\n" +
                           "(3) Press the Next> Button to accept settings and create the volume\n");
    }

    public boolean start()
    {
        return true;
    }

    public boolean end()
    {
        return true;
    }

    public int getMode()
    {
        return panel.getMode();
    }

    public StsBoxSetSubVolume[] getSetList()
    {
        return panel.getSetList();
    }
}
