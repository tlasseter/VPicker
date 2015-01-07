
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.ArbitraryLine;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;

import java.awt.*;

public class StsArbitraryLineSelect extends StsWizardStep
{
    StsArbitraryLineSelectPanel panel;
    StsHeaderPanel header;

    public StsArbitraryLineSelect(StsWizard wizard)
    {
        super(wizard);
        panel = new StsArbitraryLineSelectPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 400));
        header.setTitle("Arbitrary Line Selection");
        header.setSubtitle("Selecting Available Arbitrary Line");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#ArbitraryLine");
    }

    public StsArbitraryLine getSelectedLine()
    {
        return panel.getSelectedLine();
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

