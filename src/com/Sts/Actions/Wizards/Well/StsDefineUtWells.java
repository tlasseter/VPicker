
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Well;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.Types.*;

public class StsDefineUtWells extends StsWizardStep
{
    StsDefineUtWellsPanel panel;
    StsHeaderPanel header, info;
    String[] selectedWellnames = null;

    public StsDefineUtWells(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineUtWellsPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
//        panel.setPreferredSize(new Dimension(400, 400));
        header.setTitle("UT Well File Definition");
        header.setSubtitle("Define the origins for the selected well");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Well");
        header.setInfoText(wizardDialog,"(1)");
    }

    public void setSelectedWellnames(String[] selectedWellnames) { this.selectedWellnames = selectedWellnames; }
    public String[] getSelectedWells() { return selectedWellnames; }

    public boolean start()
    {
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        return true;
    }

    public StsPoint[] getTopHoleLocations()
    {
        return panel.getTopHoleLocations();
    }

    public byte[] getTypes()
    {
        return panel.getTypes();
    }
}

