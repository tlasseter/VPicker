
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Fracture;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

public class StsSelectFractureSet extends StsWizardStep
{
    StsSelectFractureSetPanel panel;
    StsHeaderPanel header;
    private boolean createNewFracture = false;

    public StsSelectFractureSet(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSelectFractureSetPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Select Fracture Set to Edit");
        header.setSubtitle("Selecting Available Fracture Set");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Fracture");
        header.setInfoText(wizardDialog,"(1a) Select an existing fracture set to add to.\n" +
        		"                OR                   \n" +
        		"(1b) Press the New Fracture Set Button to create a new set of fractures.\n" +
        		"(2) Press the Next>> Button");
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

