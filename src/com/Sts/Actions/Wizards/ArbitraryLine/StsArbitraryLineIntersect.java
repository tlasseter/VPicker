
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

public class StsArbitraryLineIntersect extends StsWizardStep
{
    StsArbitraryLineIntersectPanel panel;
    StsHeaderPanel header;

    public StsArbitraryLineIntersect(StsWizard wizard)
    {
        super(wizard);
        panel = new StsArbitraryLineIntersectPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(200, 400));
        header.setTitle("Arbitrary Line Definition");
        header.setSubtitle("Line Intersection Definition");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#ArbitraryLine");
    }

    public boolean start()
    {
        enableFinish();
        return true;
    }

    public boolean end()
    {
        return true;
    }

}

