
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.TimeSeries;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;

import java.awt.*;

public class StsTimeSeriesSelect extends StsWizardStep
{
    StsTimeSeriesSelectPanel panel;
    StsHeaderPanel header;

    public StsTimeSeriesSelect(StsWizard wizard)
    {
        super(wizard);
        panel = new StsTimeSeriesSelectPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 400));
        header.setTitle("Movie Selection");
        header.setSubtitle("Select Available Movie");
        header.setInfoText(wizardDialog,"(1) Select an existing animation/movie.\n" +
                           "(2) Press the Next>> Button.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#TimeSeries");
    }

    public StsMovie getSelectedMovie()
    {
        return panel.getSelectedMovie();
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

