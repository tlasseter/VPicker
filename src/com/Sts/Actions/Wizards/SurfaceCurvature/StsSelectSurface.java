
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.SurfaceCurvature;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

public class StsSelectSurface extends StsWizardStep
{
    StsSelectSurfacePanel panel;
    StsHeaderPanel header;

    public StsSelectSurface(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSelectSurfacePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Select Surface");
        header.setSubtitle("Selecting Surface for Curvature Analysis");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#SurfaceCurvature");
        header.setInfoText(wizardDialog,"(1) Select a surface to run the curvature analysis.\n" +
        		"(2) Press the Next>> Button");
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

