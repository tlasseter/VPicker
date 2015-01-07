package com.Sts.Actions.Wizards.VolumeFilter;

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

public class StsSelectVolume extends StsWizardStep
{
    StsSelectVolumePanel panel;
    StsHeaderPanel header;

    public StsSelectVolume(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSelectVolumePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 300));
        header.setTitle("Seismic Volume Filter");
        header.setSubtitle("Select Seismic Volume");
        header.setInfoText(wizardDialog,"(1) Select the volume to use as image volume.\n" +
                           "(2) Constrain the image volume by selecting a subvolume to limit the analysis.\n" +
                           "(3) Specify a name for the resulting analogue cube.\n" +
                           "(4) Press the Next>> Button.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#VolumeFilter");
    }

    public StsSeismicVolume getVolume()
    {
        return panel.getVolume();
    }

    public StsSubVolume getSubVolume()
    {
        return panel.getSubVolume();
    }
/*
    public String getName()
    {
        return panel.getCubeName();
    }
*/
    public boolean isDataFloat()
    {
        return panel.isDataFloat();
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
