
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Movie;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

import java.awt.*;

public class StsMovieVolumeDefine extends StsWizardStep
{
    StsMovieVolumeDefinePanel panel;
    StsHeaderPanel header;

    public StsMovieVolumeDefine(StsWizard wizard)
    {
        super(wizard);
        panel = new StsMovieVolumeDefinePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(200, 250));
        header.setTitle("Movie Definition");
        header.setSubtitle("PostStack3d Type Definition");
        header.setInfoText(wizardDialog, "(1) Select the type of volume animation to produce.\n" +
                           "     **** A voxel animation is currently not supported. ****\n" +
                           "     **** A point set animation will animate points by the volume they came from. ****\n" +
                           "(2) Press the Next>> Button.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Movie");
    }

    public boolean start()
    {
        return true;
    }

    public boolean end()
    {
        return true;
    }

    public byte getVolumeType()  { return panel.getVolumeType(); }

}

