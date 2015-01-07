
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

public class StsSelectWellZone extends StsWizardStep
{
    StsSelectWellZonePanel panel;
    StsHeaderPanel header;

    public StsSelectWellZone(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSelectWellZonePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Select Active Well");
        header.setSubtitle("Selecting Well and Perforation");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#ActiveWell");
        header.setInfoText(wizardDialog,"(1) Select an active well.\n" +
        		"(2) Select the active perforations. Multiple selection is allowed.\n" +
                "(3) Press Clear Current Active button to clear highlighted perfs for selected well.\n" +
                "(4) Press Clear All Active button to clear all highlighted perfs\n" +
        		"(3) Press the Finish>> Button");
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