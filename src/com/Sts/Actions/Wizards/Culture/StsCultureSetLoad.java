package com.Sts.Actions.Wizards.Culture;

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;

public class StsCultureSetLoad extends StsWizardStep implements Runnable
{
    public StsProgressPanel panel;
    private StsHeaderPanel header;
    private boolean isDone = false;
	private boolean canceled = false;

    public StsCultureSetLoad(StsWizard wizard)
    {
        super(wizard);
    }

    public void constructPanel()
    {
        panel = StsProgressPanel.constructorWithCancelButton();
        header = new StsHeaderPanel();
        super.initialize(wizard, panel, null, header);
        header.setTitle("Culture Objects Selection");
        header.setSubtitle("Process/Load Culture Object(s)");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Culture");
        header.setInfoText(wizardDialog,"(1) Once loading is complete, press the Finish Button to dismiss the screen");
    }

    public boolean start()
    {
        if(Main.isGLDebug) System.out.println("StsCultureSetLoad.start() called.");
        run();
        return true;
    }

    public void run()
    {
        if (canceled)
        {
            success = false;
            return;
        }

        try
        {
            panel.appendLine("Starting culture data loading...");
            success = ((StsCultureWizard)wizard).createCultureSet(panel);
            isDone = true;
            panel.appendLine("Culture data loaded successfully.");
            wizard.enableFinish();
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

    public boolean isDone() { return isDone; }
}
