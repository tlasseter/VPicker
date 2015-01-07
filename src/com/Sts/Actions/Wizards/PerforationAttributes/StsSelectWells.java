
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.PerforationAttributes;

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.Actions.Wizards.WizardHeaders.StsHeaderPanel;

public class StsSelectWells extends StsWizardStep
{
    StsSelectWellsPanel panel;
    StsHeaderPanel header;

    public StsSelectWells(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSelectWellsPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Perforation Correlation with Sensors");
        header.setSubtitle("Select Well");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#PerforationAttributes");
        header.setInfoText(wizardDialog,"(1) Select the well with the perforations for assignment.\n" +
        		"(2) Press the Next>> Button");
    }

    public boolean start()
    {
        panel.initialize();
        wizard.enableFinish();
        return true;
    }

    public boolean end()
    {
        return true;
    }
}