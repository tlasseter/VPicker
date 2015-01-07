
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.TimeSeries;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

import java.awt.*;

public class StsTimeSeriesDefine extends StsWizardStep
{
    StsTimeSeriesDefinePanel panel;
    StsHeaderPanel header;

    public StsTimeSeriesDefine(StsWizard wizard)
    {
        super(wizard);
        panel = new StsTimeSeriesDefinePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(200, 150));
        header.setTitle("Time Series Definition");
        header.setSubtitle("Define Time Series");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#TimeSeries");
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

