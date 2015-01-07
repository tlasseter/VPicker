
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Surfaces;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;

public class StsSurfaceLoad extends StsWizardStep implements Runnable
{
    public StsProgressPanel panel;
    private StsHeaderPanel header;

    public StsSurfaceLoad(StsWizard wizard)
    {
        super();
        panel = StsProgressPanel.constructor(5, 50);
        header = new StsHeaderPanel();
        super.initialize(wizard, panel, null, header);
        header.setTitle("Surface Selection");
        header.setSubtitle("Process/Load Surface(s)");
        header.setInfoText(wizardDialog,"(1) Once loading is complete, press the Finish Button to dismiss the screen");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Surfaces");
    }

    public boolean start()
    {
        if(Main.isGLDebug) System.out.println("StsSurfaceLoad.start() called.");
        panel.appendLine("Processing & loading selected files...");
        run();
        return true;
    }

    public void run()
    {
        try
        {
            success = ((StsSurfaceWizard)wizard).createSurfaces();
            wizard.enableFinish();
//            actionManager.endCurrentAction();
        }
        catch(Exception e)
        {
            success = false;
            StsException.outputException("StsSurfaceLoad.run() failed.", e, StsException.WARNING);
        }
    }

    public boolean end()
    {
        return true;
    }
}
