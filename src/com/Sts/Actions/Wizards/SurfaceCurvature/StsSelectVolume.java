
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.SurfaceCurvature;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

public class StsSelectVolume extends StsWizardStep
{
    StsSelectVolumePanel panel;
    StsHeaderPanel header;
    private StsVolumeCurvatureWizard wizard;

    public StsSelectVolume(StsWizard wizard)
    {
        super(wizard);
        this.wizard = (StsVolumeCurvatureWizard)wizard;
        panel = new StsSelectVolumePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Select Volume");
        header.setSubtitle("Selecting Volume for Curvature Analysis");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#VolumeCurvature");
        header.setInfoText(wizardDialog,"(1) Select a volume to run the curvature analysis.\n" +
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