
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

public class StsTypeDefine extends StsWizardStep
{
    StsTypeDefinePanel panel;
    StsHeaderPanel header;

    public StsTypeDefine(StsWizard wizard)
    {
        super(wizard);
        panel = new StsTypeDefinePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(200, 150));
        header.setTitle("Lithologic Type Definition");
        header.setSubtitle("Define Type");
        header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/LithTypes.html");
    }

    public boolean start()
    {
        disableFinish();
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        return true;
    }

    public String getTypeName()
    {
        return panel.getTypeName();
    }

    public Color getTypeColor()
    {
        return panel.getTypeColor();
    }

    public byte[] getTypeTexture()
    {
        return panel.getTypeTexture();
    }

}

