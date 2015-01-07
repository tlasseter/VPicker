
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.OSConnect;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.UI.Progress.*;

public class StsServerConnect extends StsWizardStep implements Runnable
{
    public StsProgressTextPanel panel;
    private StsHeaderPanel header;

    public StsServerConnect(StsWizard wizard)
    {
        super(wizard);
    }

    public void constructPanel()
    {
        panel = StsProgressTextPanel.constructor(5, 20);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("OpenSpirit Connect");
        header.setSubtitle("Establishing Data Server Connection");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#OSConnect");
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
            success = ((StsOSConnectWizard)wizard).establishConnection(panel);
            panel.appendLine("Connection successfully established.");
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