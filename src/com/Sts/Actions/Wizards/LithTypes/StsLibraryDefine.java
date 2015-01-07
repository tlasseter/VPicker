
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.LithTypes;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;

import java.awt.*;

public class StsLibraryDefine extends StsWizardStep
{
    StsLibraryDefinePanel panel;
    StsHeaderPanel header;

    public StsLibraryDefine(StsWizard wizard)
    {
        super(wizard);
        panel = new StsLibraryDefinePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(200, 150));
        header.setTitle("Lithologic Library Definition");
        header.setSubtitle("Define Library");
        header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/LithLibrary.html");
    }

    public boolean start()
    {
        disableFinish();
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        if(panel.getLibraryName() != null)
            StsTypeLibrary.getCreate(panel.getLibraryName());
        return true;
    }

    public String getLibraryName()
    {
        return panel.getLibraryName();
    }
}

