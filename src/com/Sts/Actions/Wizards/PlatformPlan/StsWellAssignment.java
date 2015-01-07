
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.PlatformPlan;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

import java.awt.*;

public class StsWellAssignment extends StsWizardStep
{
    StsWellAssignmentPanel panel;
    StsHeaderPanel header, info;

    public StsWellAssignment(StsWizard wizard)
    {
        super(wizard);
        panel = new StsWellAssignmentPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 700));
        header.setTitle("Well Assignment");
        header.setSubtitle("Assign existing wells to slots");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#PlatformPlan");
        header.setInfoText(wizardDialog,"(1) ....\n");
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

