
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.ActiveWell;

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.Actions.Wizards.ProximityAnalysis.StsSelectSensorsPanel;
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
        header.setTitle("Select Active Well");
        header.setSubtitle("Select Active Wells");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#ActiveWell");
        header.setInfoText(wizardDialog,"(1) Press the Deactivate All Button to clear all highlighted wells\n" +
        		"(2) Select one or more wells to highlight.\n" +
        		"(3) Press the Next>> Button");
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