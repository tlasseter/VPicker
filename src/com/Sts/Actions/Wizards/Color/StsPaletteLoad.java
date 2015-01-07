
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Color;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.UI.Progress.*;

public class StsPaletteLoad extends StsWizardStep implements Runnable
{
    public StsProgressTextPanel panel;
    private StsHeaderPanel header;

    private double progressValue = 0.0;
    private boolean isDone = false;

    public StsPaletteLoad(StsWizard wizard)
    {
        super(wizard);
    }

    public void constructPanel()
    {
        panel = StsProgressTextPanel.constructor(5, 20);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Palette Selection");
        header.setSubtitle("Load Palette(s)");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Color");
    }

    public boolean start()
    {
        run();
        return true;
    }

    public void run()
    {
        try
        {
            panel.appendLine("Starting palette loading...");
            success = ((StsPaletteWizard)wizard).createSpectrums(panel);
            progressValue = 100.0f;
            isDone = true;
            panel.appendLine("Palettes loaded successfully.");
            if(!success) return;
        }
        catch(Exception e)
        {
            success = false;
        }
    }

    public boolean end()
    {
        return true;
    }
}