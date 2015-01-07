
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

public class StsDefineWell extends StsWizardStep
{
    StsDefineWellPanel panel;
    StsHeaderPanel header;

    public StsDefineWell(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineWellPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Create a Real-time Monitor");
        header.setSubtitle("Define a New Well to Monitor");
        header.setInfoText(wizardDialog, "(1) Supply the well name.\n" +
                           "(8) Press the Next>> Button.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Monitor");
    }

    public boolean start()
    {
        wizard.enableFinish();
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        return true;
    }

}