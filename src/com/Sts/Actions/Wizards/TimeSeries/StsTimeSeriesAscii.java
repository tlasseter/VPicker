
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

public class StsTimeSeriesAscii extends StsWizardStep
{
    StsTimeSeriesAsciiPanel panel;
    StsHeaderPanel header;

    public StsTimeSeriesAscii(StsWizard wizard)
    {
        super(wizard);
        panel = new StsTimeSeriesAsciiPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 400));
        header.setTitle("Time Series Import");
        header.setSubtitle("Import Ascii Time Series File");
        header.setInfoText(wizardDialog, "(1) Press the Select button to find and import Ascii file.\n" +
                           "(2) Press the Next>> Button.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#TimeSeries");
    }

    public boolean start()
    {
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        panel.determineTimeType();
        return true;
    }

}

