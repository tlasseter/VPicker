
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.VolumeStimulated;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

public class StsSelectSubVolumes extends StsWizardStep
{
    StsSelectSubVolumesPanel panel;
    StsHeaderPanel header;

    public StsSelectSubVolumes(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSelectSubVolumesPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Select SubVolumes");
        header.setSubtitle("Select sub-volumes to limit stimulated volume");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#VolumeStimulated");
        header.setInfoText(wizardDialog,"(1) Select sub-volumes from list.\n" +
        		"(2) Press the Next>> Button");
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

