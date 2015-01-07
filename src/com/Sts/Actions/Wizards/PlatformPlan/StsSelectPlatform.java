
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.PlatformPlan;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;

import java.awt.*;

public class StsSelectPlatform extends StsWizardStep
{
    StsSelectPlatformPanel panel;
    StsHeaderPanel header, info;

    public StsSelectPlatform(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSelectPlatformPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 700));
        header.setTitle("Platform Selection");
        header.setSubtitle("Selecting Available Platform");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#PlatformPlan");
        header.setInfoText(wizardDialog,"(1) Select an existing platform to edit. \n" +
                           "                   OR                    \n" +
                           "(2) Select New Platform to create a new platform\n");
    }

    public StsPlatform getSelectedPlatform()
    {
        return panel.getSelectedPlatform();
    }

    public boolean start()
    {
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        panel.getSelectedPlatform();
        return true;
    }
}

