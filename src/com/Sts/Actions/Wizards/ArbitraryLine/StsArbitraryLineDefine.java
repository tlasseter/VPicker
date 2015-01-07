
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.ArbitraryLine;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

import java.awt.*;

public class StsArbitraryLineDefine extends StsWizardStep
{
    StsArbitraryLineDefinePanel panel;
    StsHeaderPanel header;

    public StsArbitraryLineDefine(StsWizard wizard)
    {
        super(wizard);
        panel = new StsArbitraryLineDefinePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(200, 400));
        header.setTitle("Arbitrary Line Definition");
        header.setSubtitle("Define Arbitrary Line");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#ArbitraryLine");
    }

    public boolean start()
    {
        disableFinish();
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        panel.activateRotation(0.0f);
        return true;
    }
}

