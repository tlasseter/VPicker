
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.MakeMovie;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.UI.Progress.*;

public class StsMakeMovieCreate extends StsWizardStep implements Runnable
{
    public StsMakeMovieWizard wizard;
    public StsProgressPanel panel;
    private StsHeaderPanel header;

    private boolean isDone = false;
	private boolean canceled = false;

    public StsMakeMovieCreate(StsWizard wizard)
    {
        super(wizard);
        this.wizard = (StsMakeMovieWizard)wizard;
    }

    public void constructPanel()
    {
        panel = StsProgressPanel.constructorWithCancelButton();
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Quicktime Movie Definition");
        header.setSubtitle("Create Movie");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#MakeMovie");                                
    }

    public boolean start()
    {
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
            disableFinish();

            // Setup the Movie Object
            panel.setValue(0.0f);
            panel.setDescription("Configuring Movie: " + wizard.getOutputMovieName());
            if(wizard.createMovie(panel))
                panel.appendLine("Successfully created movie.");
            else
                panel.appendLine("Failed to created movie.");

            panel.setValue(100.0f);
            enableFinish();
            disableCancel();
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

    public boolean isDone() { return isDone; }
}
