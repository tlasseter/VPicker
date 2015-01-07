
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Pointset;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.UI.*;

public class StsPointsetLoad extends StsWizardStep implements Runnable
{
    public StsStatusPanel panel;
    private StsHeaderPanel header;

    public StsPointsetLoad(StsWizard wizard)
    {
        super(wizard, new StsStatusPanel(), null, new StsHeaderPanel());
        panel = (StsStatusPanel) getContainer();

        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Pointset Selection");
        header.setSubtitle("Load Pointset(s)");
        header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/PointsetLoad.html");
    }

    public boolean start()
    {
        panel.setTitle("Loading selected files...");

        run();
        return true;
    }

    public void run()
    {
        try
        {
            success = ((StsPointsetWizard)wizard).createPointsets();
            actionManager.endCurrentAction();
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
