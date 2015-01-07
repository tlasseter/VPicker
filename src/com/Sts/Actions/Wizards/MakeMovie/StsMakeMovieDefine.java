
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.MakeMovie;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

import java.awt.*;

public class StsMakeMovieDefine extends StsWizardStep
{
    StsMakeMovieDefinePanel panel;
    StsHeaderPanel header;

    public StsMakeMovieDefine(StsWizard wizard)
    {
        super(wizard);
        panel = new StsMakeMovieDefinePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(300, 400));
        header.setTitle("Quicktime Movie Definition");
        header.setSubtitle("Capture and Define Movie");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#MakeMovie");                                        
        header.setInfoText(wizardDialog, "(1)   \n" +
                           "(2)   \n" +
                           "(3) Press the Next>> Button.\n");
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

