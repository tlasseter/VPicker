
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.CrossPlot;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;

import java.awt.*;

public class StsSelectCrossplot extends StsWizardStep
{
    StsSelectCrossplotPanel panel;
    StsHeaderPanel header;
    private boolean createNewCrossplot = false;

    public StsSelectCrossplot(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSelectCrossplotPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 300));
        header.setTitle("Cross Plot Selection");
        header.setSubtitle("Selecting Available Cross Plots");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Crossplot");
    }

    public StsCrossplot getSelectedCrossplot()
    {
        return panel.getSelectedCrossplot();
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

