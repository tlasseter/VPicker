
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.SubVolumeMgmt;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

import java.awt.*;

public class StsSubVolumeMgmtActivate extends StsWizardStep
{
    StsSubVolumeMgmtActivatePanel panel;
    StsHeaderPanel header;

    public StsSubVolumeMgmtActivate(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSubVolumeMgmtActivatePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(300, 400));
        header.setTitle("Sub-PostStack3d Management");
        header.setSubtitle("Activate Sub-Volumes");
        header.setInfoText(wizardDialog, "(1)   \n" +
                           "(2)   \n" +
                           "(3) Press the Next>> Button.\n");
        header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/MakeMovie.html");
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

