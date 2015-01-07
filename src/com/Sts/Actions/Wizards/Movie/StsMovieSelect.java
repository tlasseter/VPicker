
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Movie;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;

import java.awt.*;

public class StsMovieSelect extends StsWizardStep
{
    StsMovieSelectPanel panel;
    StsHeaderPanel header;
    private boolean createNewCrossplot = false;

    public StsMovieSelect(StsWizard wizard)
    {
        super(wizard);
        panel = new StsMovieSelectPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 300));
        header.setTitle("Movie Selection");
        header.setSubtitle("Selecting Available Movies");
        header.setInfoText(wizardDialog,"(1) Select an existing animation/movie or press the New Movie button.\n" +
                           "     **** The New Movie button will automatically proceed to the next step. *****\n" +
                           "     **** If selecting an existing animation/movie the Next>> button must be pressed.*****\n" +
                           "(2) Press the Next>> Button.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Movie");
    }

    public StsMovie getSelectedMovie()
    {
        return panel.getSelectedMovie();
    }

    public boolean start()
    {
        return true;
    }

    public boolean end()
    {
        return true;
    }
}

