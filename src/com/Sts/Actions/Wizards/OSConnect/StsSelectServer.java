
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.OSConnect;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

import java.awt.*;

public class StsSelectServer extends StsWizardStep
{
	StsSelectServerPanel panel;
    StsHeaderPanel header;

    public StsSelectServer(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSelectServerPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 300));
        header.setTitle("Server Selection");
        header.setSubtitle("Select an OpenSpirit Server");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#OSConnect");
        header.setInfoText(wizardDialog,"(1) Select one of the available OpenSpirit Servers\n" +
                           "(2) Once desired server has been selected press the Next>> Button");
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

