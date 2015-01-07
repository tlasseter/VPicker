
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.OSWell;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

public class StsOSWellSelect extends StsWizardStep
{
    StsOSWellSelectPanel panel;
    StsHeaderPanel header, info;

    public StsOSWellSelect(StsWizard wizard)
    {
        super(wizard);
        panel = new StsOSWellSelectPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("OpenSpirit Well Selection");
        header.setSubtitle("Select from Available Wells");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#OSWell");
        header.setInfoText(wizardDialog,"(1) Select the desired wells and place them in the right list using the provided controls\n " +
                           "(2) Once well selections are complete, press the Next>> Button.");
    }

    public StsOSWell[] getSelectedWells() {  return panel.getSelectedWells(); }
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

