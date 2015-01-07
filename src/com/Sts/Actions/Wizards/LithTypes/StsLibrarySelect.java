
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.LithTypes;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

import java.awt.*;

public class StsLibrarySelect extends StsWizardStep
{
    StsLibrarySelectPanel panel;
    StsHeaderPanel header;
    private boolean createNewCrossplot = false;

    public StsLibrarySelect(StsWizard wizard)
    {
        super(wizard);
        panel = new StsLibrarySelectPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 300));
        header.setTitle("Library Selection");
        header.setSubtitle("Selecting Available Library");
        header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/LithTypes.html");
    }

    public boolean start()
    {
        return true;
    }

    public boolean end()
    {
        ((StsLithTypesWizard)wizard).setSelectedLibrary(panel.getSelectedLibrary());
        return true;
    }
}

