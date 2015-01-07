
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.LogToVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;

public class StsSelectWells extends StsWizardStep
{
    StsSelectWellsPanel panel;
    StsHeaderPanel header;
    StsWell[] wells = null;
    int[] selected = null;

    public StsSelectWells(StsWizard wizard)
    {
        super(wizard);
        panel = new StsSelectWellsPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Create Volume from Log");
        header.setSubtitle("Select Wells");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#LogToVolume");                        
        header.setInfoText(wizardDialog,"(1) Select all wells to use in volume construction.\n" +
                                  "   ***** the log will be selected from the selected wells ***** \n" +
                                  "(2) Press the Next>> Button.");
    }

    public boolean start()
    {
        panel.initialize();
        return true;
    }

    public boolean end() { return getSelectedWells() != null; }

    public Object[] getSelectedWells()
    {
        return panel.getSelectedWells();
    }
}