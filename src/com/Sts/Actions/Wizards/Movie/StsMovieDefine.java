
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

public class StsMovieDefine extends StsWizardStep
{
    StsMovieDefinePanel panel;
    StsHeaderPanel header;

    public StsMovieDefine(StsWizard wizard)
    {
        super(wizard);
        panel = new StsMovieDefinePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(300, 400));
        header.setTitle("Movie Definition");
        header.setSubtitle("Define Movie");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Movie");
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
    public float[] getRange()
    {
        return panel.getRange();
    }
    public float getIncrement()
    {
        return panel.getIncrement();
    }
    public int getDelay()
    {
        return panel.getDelay();
    }
    public boolean getLoop()
    {
        return panel.getLoop();
    }
    public boolean getCycleVolumes()
    {
        return panel.getCycleVolumes();
    }
    public boolean getElevation()
    {
        return panel.getElevation();
    }
    public boolean getAzimuth()
    {
        return panel.getAzimuth();
    }
    public int getAzimuthStart()
    {
        return panel.getAzimuthStart();
    }
    public int getAzimuthIncrement()
    {
        return panel.getAzimuthIncrement();
    }
    public int getElevationStart()
    {
        return panel.getElevationStart();
    }
    public int getElevationIncrement()
    {
        return panel.getElevationIncrement();
    }

}

