
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

public class StsMovieDataDefine extends StsWizardStep
{
    StsMovieDataDefinePanel panel;
    StsHeaderPanel header;

    public StsMovieDataDefine(StsWizard wizard)
    {
        super(wizard);
        panel = new StsMovieDataDefinePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(200, 250));
        header.setTitle("Movie Definition");
        header.setSubtitle("Line Type Definition");
        header.setInfoText(wizardDialog,"(1) Specify the range and increment.\n" +
                           "(2) Specify the elevation (vertical rotation) animation parameters.\n" +
                           "     **** Check elevation check box to rotate vertically. ****\n" +
                           "     **** Specify the starting elevation of user in degrees. ****\n" +
                           "     **** Specify the increment in degrees to change elevation each frame. ****\n" +
                           "(3) Specify the azimuth (horizontal rotation) animation parameters.\n" +
                           "     **** Check azimuth check box to rotate horizontally. ****\n" +
                           "     **** Specify the starting azimuth of user in degrees. ****\n" +
                           "     **** Specify the increment in degrees to change azimuth each frame. ****\n" +
                           "(4) Specify the frame rate in milliseconds, or selected full speed check box.\n" +
                           "(5) Specify whether to loop after end of animation is reached.\n" +
                           "(6) Press the Next>> Button to build animation and activate movie toolbar.\n");

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

    public int getDirection()
    {
        return panel.getDirection();
    }

}

