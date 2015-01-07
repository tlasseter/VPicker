package com.Sts.Actions.Wizards.MicroseismicPreStack;

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

public class StsMicroPreStackLoad extends StsWizardStep implements Runnable
{
    public StsProgressPanel panel;
    private StsHeaderPanel header;
    private boolean isDone = false;
	private boolean canceled = false;

    public StsMicroPreStackLoad(StsWizard wizard)
    {
        super(wizard);
    }

    public void constructPanel()
    {
        panel = StsProgressPanel.constructorWithCancelButton(10, 50);
        header = new StsHeaderPanel();
        setPanels(panel, header);    	
        header.setTitle("Microseismic PreStack Load");
        header.setSubtitle("Process/Load Selected Files");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#MicroseismicPreStack");
        header.setInfoText(wizardDialog,"(1) Once loading is complete, press the Finish Button to dismiss the screen");
    }

    public boolean start()
    {
        if(Main.isGLDebug) System.out.println("StsMicroPreStackLoad.start() called.");
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
            panel.appendLine("Starting microseismic prestack data loading...");
            success = ((StsMicroPreStackWizard)wizard).createMicroseismicSet(panel);
            isDone = true;
            panel.appendLine("Microseismic prestack data loaded successfully.");
            wizard.enableFinish();
        }
        catch(Exception e)
        {
            success = false;
            StsException.outputException("StsMicroPreStackLoad.run() failed.", e, StsException.WARNING);
        }
    }

    public boolean end()
    {
        return true;
    }

    public boolean isDone() { return isDone; }
}
