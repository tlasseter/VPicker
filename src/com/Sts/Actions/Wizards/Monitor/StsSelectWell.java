
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Monitor;

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.Actions.Wizards.WizardHeaders.StsHeaderPanel;

public class StsSelectWell extends StsWizardStep
{
    StsSelectWellPanel panel;
    StsHeaderPanel header;

    public StsSelectWell(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSelectWellPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Create a Real-time Monitor");
        header.setSubtitle("Select Well to Monitor");
        header.setInfoText(wizard.dialog, "(1) Select one well from the list.\n" +
        		"                 OR\n" +
        		"(2) Press the New Drilling Well... Button to create a new well.\n" +
            "(3) Press the Next Button.");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Monitor");
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